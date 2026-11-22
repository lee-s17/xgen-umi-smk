# xGen UMI duplex pipeline (SEQC2 ctDNA, IDT assay)
#
# FASTQ -> fastp -> FastqToBam (UMI) -> bwa | ZipperBams -> GroupReadsByUmi
#   -> Call{Duplex,Molecular}ConsensusReads -> bwa | ZipperBams
#   -> FilterConsensusReads (one pass per family-size level) -> ClipBam
#   -> VarDict / Mutect2 -> bcftools filter -> benchmark vs truth

import os
import re
import pandas as pd

configfile: "config/config.yaml"

OUT = config["outdir"]
REF = config["ref"]
REF_DICT = re.sub(r"\.(fa|fasta|fna)(\.gz)?$", ".dict", REF)
CALLER = config["caller"]
MODE = config["consensus"]["mode"]
LEVELS = config["filter_levels"]

samples = (
    pd.read_csv(config["samples"], sep="\t", comment="#", dtype=str)
    .fillna("")
    .set_index("sample", drop=False)
)
if "role" not in samples.columns:
    samples["role"] = "test"
samples["role"] = samples["role"].replace("", "test")
SAMPLES = list(samples.index)
TEST_SAMPLES = [s for s in SAMPLES if samples.loc[s, "role"] != "normal"]
NORMALS = [s for s in SAMPLES if samples.loc[s, "role"] == "normal"]
BENCH_SAMPLES = [s for s in TEST_SAMPLES if samples.loc[s, "truth"]]
MAX_TLEN = int(config.get("fragment", {}).get("max_tlen", 0) or 0)

wildcard_constraints:
    sample="[^/.]+",
    level="[^/.]+",
    truth="[^/.]+",


def fgbio(mem_gb):
    # bioconda fgbio wrapper passes JVM flags through; level-1 compression for intermediates
    return lambda w, resources: (
        f"fgbio -Xmx{mem_gb}g -Djava.io.tmpdir={resources.tmpdir} --compression 1 --async-io true"
    )


def raw_fastq(w):
    s = samples.loc[w.sample]
    if s["fq1"]:
        return {"r1": s["fq1"], "r2": s["fq2"]}
    return {"r1": f"{OUT}/raw/{w.sample}_R1.fastq.gz", "r2": f"{OUT}/raw/{w.sample}_R2.fastq.gz"}


def min_reads(level):
    m = LEVELS[level]["min_reads"]
    # simplex consensus takes a single integer
    return m if MODE == "duplex" else m.split()[0]


rule all:
    input:
        expand(f"{OUT}/consensus/{{sample}}.{{level}}.bam", sample=SAMPLES, level=LEVELS),
        expand(f"{OUT}/vcf/{{sample}}.{{level}}.vcf.gz", sample=TEST_SAMPLES, level=LEVELS),
        expand(f"{OUT}/qc/fragments/{{sample}}.{{level}}.stats.txt", sample=SAMPLES, level=LEVELS),
        expand(f"{OUT}/informed/{{sample}}.{{level}}.summary.tsv", sample=BENCH_SAMPLES, level=LEVELS),
        expand(f"{OUT}/metrics/{{sample}}.duplex_yield_metrics.txt", sample=SAMPLES)
        if MODE == "duplex" else [],
        f"{OUT}/qc/multiqc_report.html",
        [f"{OUT}/benchmark/summary.tsv"] if BENCH_SAMPLES else [],


# Reference and panel

rule ref_fai:
    input: ancient(REF)
    output: REF + ".fai"
    conda: "envs/core.yaml"
    shell: "samtools faidx {input}"


rule ref_dict:
    input: ancient(REF)
    output: REF_DICT
    conda: "envs/core.yaml"
    shell: "samtools dict {input} -o {output}"


# The index is not an input of the alignment rules, so a moved/rebuilt index never
# re-triggers alignment. Build it once with: snakemake bwa_index_all
rule bwa_index_all:
    input: REF + ".bwt.2bit.64"


rule bwa_index:
    # bwa-mem2 index: ~5x larger than bwa's and needs tens of GB RAM to build; built once per reference
    input: ancient(REF)
    output: REF + ".bwt.2bit.64"
    resources: mem_mb=70000
    log: f"{OUT}/logs/bwa_index.log"
    conda: "envs/core.yaml"
    shell: "bwa-mem2 index {input} &> {log}"


rule prep_bed:
    # sort/merge, add a name column (VarDict -g 4) and build a Picard interval_list
    input:
        bed=config["panel_bed"],
        dict=ancient(REF_DICT),
    output:
        bed=f"{OUT}/ref/panel.bed",
        ilist=f"{OUT}/ref/panel.interval_list",
    conda: "envs/core.yaml"
    shell:
        """
        grep -vE '^(track|browser|#)' {input.bed} | cut -f1-3 \
          | sort -k1,1 -k2,2n | bedtools merge -i - \
          | awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$1":"$2+1"-"$3}}' > {output.bed}
        cat {input.dict} > {output.ilist}
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$2+1,$3,"+",$4}}' {output.bed} >> {output.ilist}
        """


# Download (only for samples with an SRR and no local FASTQ)

rule fetch_sra:
    # temp(): deleted after fastqc/fastp, and never re-downloaded when later outputs exist
    output:
        r1=temp(f"{OUT}/raw/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{OUT}/raw/{{sample}}_R2.fastq.gz"),
    params:
        # several runs per library are allowed, separated by ';' (concatenated in order)
        srr=lambda w: " ".join(samples.loc[w.sample, "srr"].split(";")),
        tmp=lambda w: f"{OUT}/raw/{w.sample}_sra",
    threads: 6
    log: f"{OUT}/logs/fetch_sra/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        (set -euo pipefail
        mkdir -p {params.tmp}
        : > {output.r1}; : > {output.r2}
        for r in {params.srr}; do
            prefetch $r -O {params.tmp} --max-size 200G
            (cd {params.tmp} && fasterq-dump $r --split-files --skip-technical -e {threads})
            pigz -p {threads} -c {params.tmp}/${{r}}_1.fastq >> {output.r1}
            pigz -p {threads} -c {params.tmp}/${{r}}_2.fastq >> {output.r2}
            rm -f {params.tmp}/${{r}}_*.fastq
        done
        rm -rf {params.tmp}) &> {log}
        """


# QC and trimming

rule fastqc:
    input: unpack(raw_fastq)
    output: directory(f"{OUT}/qc/fastqc/{{sample}}")
    threads: 2
    log: f"{OUT}/logs/fastqc/{{sample}}.log"
    conda: "envs/core.yaml"
    shell: "mkdir -p {output} && fastqc --quiet --threads {threads} --outdir {output} {input.r1} {input.r2} &> {log}"


rule fastp:
    # 3'-only trimming (adapter, quality tail, poly-G) so the 5' UMI+skip stays in place.
    # No base correction (-c): it would rewrite bases before consensus.
    input: unpack(raw_fastq)
    output:
        r1=temp(f"{OUT}/trimmed/{{sample}}_R1.fastq.gz"),
        r2=temp(f"{OUT}/trimmed/{{sample}}_R2.fastq.gz"),
        json=f"{OUT}/qc/fastp/{{sample}}.fastp.json",
        html=f"{OUT}/qc/fastp/{{sample}}.fastp.html",
    params: extra=config["fastp"]["extra"]
    threads: 8
    log: f"{OUT}/logs/fastp/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
          -w {threads} -z 1 {params.extra} \
          -j {output.json} -h {output.html} &> {log}
        """


# UMI extraction straight from FASTQ (replaces FastqToSam + ExtractUmisFromBam)

rule fastq_to_ubam:
    input:
        r1=f"{OUT}/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/trimmed/{{sample}}_R2.fastq.gz",
    output: temp(f"{OUT}/bam/{{sample}}.unmapped.bam")
    params:
        rs=config["read_structure"],
        fg=fgbio(4),
    resources: mem_mb=6000
    log: f"{OUT}/logs/fastq_to_ubam/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        {params.fg} FastqToBam -i {input.r1} {input.r2} -r {params.rs} \
          --umi-tag RX --sample {wildcards.sample} --library {wildcards.sample} \
          --read-group-id {wildcards.sample} -o {output} &> {log}
        """


# Raw alignment: streamed, no intermediate FASTQ/BAM, no MergeBamAlignment

rule align_raw:
    input:
        ubam=f"{OUT}/bam/{{sample}}.unmapped.bam",
        dict=ancient(REF_DICT),
        fai=ancient(REF + ".fai"),
    output: temp(f"{OUT}/bam/{{sample}}.mapped.tc.bam")
    params:
        fg=fgbio(4),
        sort_threads=4,
    threads: 16
    resources: mem_mb=24000
    log: f"{OUT}/logs/align_raw/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        test -e {REF}.bwt.2bit.64 || {{ echo "missing bwa-mem2 index: run snakemake bwa_index_all or bwa-mem2 index {REF}" >&2; exit 1; }}
        (samtools fastq {input.ubam} \
          | bwa-mem2 mem -t {threads} -p -K 150000000 -Y {REF} - \
          | {params.fg} ZipperBams --unmapped {input.ubam} --ref {REF} --output /dev/stdout \
          | samtools sort --template-coordinate -@ {params.sort_threads} -m 2G \
              -T {resources.tmpdir}/{wildcards.sample}.tc -o {output}) 2> {log}
        """


rule group_umis:
    input: f"{OUT}/bam/{{sample}}.mapped.tc.bam"
    output:
        bam=temp(f"{OUT}/bam/{{sample}}.grouped.bam"),
        hist=f"{OUT}/metrics/{{sample}}.family_size_histogram.txt",
    params:
        fg=fgbio(8),
        g=config["group"],
    resources: mem_mb=10000
    log: f"{OUT}/logs/group_umis/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        {params.fg} GroupReadsByUmi -i {input} -o {output.bam} -f {output.hist} \
          --strategy {params.g[strategy]} --edits {params.g[edits]} \
          --min-map-q {params.g[min_map_q]} -t RX -T MI &> {log}
        """


rule duplex_metrics:
    # family-size distributions + duplex yield downsampling curve (plateau analysis)
    input:
        bam=f"{OUT}/bam/{{sample}}.grouped.bam",
        ilist=f"{OUT}/ref/panel.interval_list",
    output:
        f"{OUT}/metrics/{{sample}}.duplex_yield_metrics.txt",
        f"{OUT}/metrics/{{sample}}.family_sizes.txt",
        f"{OUT}/metrics/{{sample}}.duplex_family_sizes.txt",
    params:
        fg=fgbio(8),
        prefix=f"{OUT}/metrics/{{sample}}",
    resources: mem_mb=12000
    log: f"{OUT}/logs/duplex_metrics/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        {params.fg} CollectDuplexSeqMetrics -i {input.bam} -o {params.prefix} \
          -l {input.ilist} --duplex-umi-counts true &> {log}
        """


rule call_consensus:
    # deliberately permissive (min-reads 1 [1 0]); stringency is applied later by
    # FilterConsensusReads so several family-size levels reuse one consensus pass
    input: f"{OUT}/bam/{{sample}}.grouped.bam"
    output: temp(f"{OUT}/bam/{{sample}}.consensus.unmapped.bam")
    params:
        fg=fgbio(8),
        c=config["consensus"],
        tool="CallDuplexConsensusReads" if MODE == "duplex" else "CallMolecularConsensusReads",
        minr="1 1 0" if MODE == "duplex" else "1",
    threads: 8
    resources: mem_mb=12000
    log: f"{OUT}/logs/call_consensus/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        {params.fg} {params.tool} -i {input} -o {output} \
          --error-rate-pre-umi {params.c[error_rate_pre_umi]} \
          --error-rate-post-umi {params.c[error_rate_post_umi]} \
          --min-input-base-quality {params.c[min_input_base_quality]} \
          --min-reads {params.minr} --threads {threads} &> {log}
        """


rule align_consensus:
    # output stays query-grouped: FilterConsensusReads/ClipBam iterate by template
    input:
        ubam=f"{OUT}/bam/{{sample}}.consensus.unmapped.bam",
        dict=ancient(REF_DICT),
    output: f"{OUT}/bam/{{sample}}.consensus.mapped.bam"
    params: fg=fgbio(4)
    threads: 16
    resources: mem_mb=16000
    log: f"{OUT}/logs/align_consensus/{{sample}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        test -e {REF}.bwt.2bit.64 || {{ echo "missing bwa-mem2 index: run snakemake bwa_index_all or bwa-mem2 index {REF}" >&2; exit 1; }}
        (samtools fastq {input.ubam} \
          | bwa-mem2 mem -t {threads} -p -K 150000000 -Y {REF} - \
          | {params.fg} ZipperBams --unmapped {input.ubam} --ref {REF} \
              --tags-to-reverse Consensus --tags-to-revcomp Consensus \
              --output {output}) 2> {log}
        """


rule filter_clip:
    # one cheap pass per family-size level; overlap clipping avoids double-counting
    # a molecule's R1/R2 overlap in VAF
    # FilterConsensusReads reads the ZipperBams output as-is (query-grouped; a samtools
    # name sort breaks its R1/R2 pairing). fgbio 4.x ClipBam insists on SO:queryname, so
    # the filtered stream is sorted with fgbio's own SortBam before clipping.
    input: f"{OUT}/bam/{{sample}}.consensus.mapped.bam"
    output:
        bam=f"{OUT}/consensus/{{sample}}.{{level}}.bam",
        bai=f"{OUT}/consensus/{{sample}}.{{level}}.bam.bai",
        clip=f"{OUT}/metrics/{{sample}}.{{level}}.clip_metrics.txt",
    params:
        fg=fgbio(4),
        f=config["filter"],
        minr=lambda w: min_reads(w.level),
        # optional ctDNA-enrichment filter (INVAR): keep templates <= max_tlen; unmapped/unpaired keep tlen 0
        frag=(f"samtools view -u -e 'tlen==0 || abs(tlen)<={MAX_TLEN}' -" if MAX_TLEN else "cat"),
    threads: 4
    resources: mem_mb=22000   # 3 JVMs (4g heap each, ~5.5 GB real) + samtools sort
    log: f"{OUT}/logs/filter_clip/{{sample}}.{{level}}.log"
    conda: "envs/core.yaml"
    shell:
        """
        ({params.fg} FilterConsensusReads -i {input} -o /dev/stdout --ref {REF} \
          --min-reads {params.minr} \
          --min-base-quality {params.f[min_base_quality]} \
          --max-base-error-rate {params.f[max_base_error_rate]} \
          --max-read-error-rate {params.f[max_read_error_rate]} \
          --max-no-calls {params.f[max_no_calls]} \
          --reverse-per-base-tags true \
        | {params.fg} SortBam -i /dev/stdin -o /dev/stdout -s Queryname --max-records-in-ram 1000000 \
        | {params.fg} ClipBam -i /dev/stdin -o /dev/stdout --ref {REF} \
          --clip-overlapping-reads true -m {output.clip} \
        | {params.frag} \
        | samtools sort -@ {threads} -m 1G -T {resources.tmpdir}/{wildcards.sample}.{wildcards.level} \
          -o {output.bam}
        samtools index {output.bam}) 2> {log}
        """


rule mosdepth:
    input:
        bam=f"{OUT}/consensus/{{sample}}.{{level}}.bam",
        bed=f"{OUT}/ref/panel.bed",
    output: f"{OUT}/qc/mosdepth/{{sample}}.{{level}}.mosdepth.summary.txt"
    params: prefix=f"{OUT}/qc/mosdepth/{{sample}}.{{level}}"
    threads: 4
    log: f"{OUT}/logs/mosdepth/{{sample}}.{{level}}.log"
    conda: "envs/core.yaml"
    shell: "mosdepth -n --fast-mode -t {threads} -b {input.bed} {params.prefix} {input.bam} &> {log}"


rule fragment_stats:
    # insert-size distribution of consensus templates (MultiQC picks up samtools stats)
    input: f"{OUT}/consensus/{{sample}}.{{level}}.bam"
    output: f"{OUT}/qc/fragments/{{sample}}.{{level}}.stats.txt"
    threads: 2
    conda: "envs/core.yaml"
    shell: "samtools stats -@ {threads} {input} > {output}"


# Variant calling

if CALLER == "vardict":

    rule vardict:
        input:
            bam=f"{OUT}/consensus/{{sample}}.{{level}}.bam",
            bed=f"{OUT}/ref/panel.bed",
        output: temp(f"{OUT}/vcf/{{sample}}.{{level}}.raw.vcf")
        params: v=config["vardict"]
        threads: 8
        resources: mem_mb=16000
        log: f"{OUT}/logs/vardict/{{sample}}.{{level}}.log"
        conda: "envs/callers.yaml"
        shell:
            """
            (set -euo pipefail
            export JAVA_OPTS="-Xmx12g"
            raw={output}.vardict.txt
            vardict-java -G {REF} -N {wildcards.sample} -b {input.bam} \
              -f {params.v[min_af]} -q {params.v[min_bq]} -Q {params.v[min_mq]} \
              -r 1 -th {threads} -c 1 -S 2 -E 3 -g 4 {params.v[extra]} {input.bed} > $raw
            # teststrandbias.R stops with "no lines available" on empty input, so only run it
            # when VarDict reported something; var2vcf still writes a header-only VCF
            if [ -s $raw ]; then
                teststrandbias.R < $raw | var2vcf_valid.pl -N {wildcards.sample} -E -f {params.v[min_af]} > {output}
            else
                var2vcf_valid.pl -N {wildcards.sample} -E -f {params.v[min_af]} < /dev/null > {output}
            fi
            rm -f $raw
            grep -q '^##fileformat=VCF' {output}) 2> {log}
            """

elif CALLER == "mutect2":

    rule mutect2:
        input:
            bam=f"{OUT}/consensus/{{sample}}.{{level}}.bam",
            ilist=f"{OUT}/ref/panel.interval_list",
        output: temp(f"{OUT}/vcf/{{sample}}.{{level}}.raw.vcf")
        params:
            m=config["mutect2"],
            unf=f"{OUT}/vcf/{{sample}}.{{level}}.m2.unfiltered.vcf.gz",
        threads: 4
        resources: mem_mb=12000
        log: f"{OUT}/logs/mutect2/{{sample}}.{{level}}.log"
        conda: "envs/callers.yaml"
        shell:
            """
            (gatk --java-options "-Xmx8g" Mutect2 -R {REF} -I {input.bam} -L {input.ilist} \
              --max-reads-per-alignment-start 0 --callable-depth {params.m[callable_depth]} \
              --minimum-allele-fraction 0 --max-mnp-distance 0 \
              --native-pair-hmm-threads {threads} {params.m[extra]} -O {params.unf}
            gatk --java-options "-Xmx4g" FilterMutectCalls -R {REF} -V {params.unf} \
              --min-allele-fraction {params.m[min_af]} -O {output}) 2> {log}
            """

else:
    raise ValueError(f"caller must be vardict or mutect2, got {CALLER}")


rule filter_vcf:
    # add ##contig lines (VarDict writes none, so BCF streams fail), normalise (-a splits
    # MNVs into SNVs so they match the per-base normal pileup and truth; split
    # multi-allelics, left-align), drop caller filters that are not meaningful on
    # consensus reads, then apply per-level alt-support threshold
    input: f"{OUT}/vcf/{{sample}}.{{level}}.raw.vcf"
    output:
        vcf=temp(f"{OUT}/vcf/{{sample}}.{{level}}.filt.vcf.gz"),
        tbi=temp(f"{OUT}/vcf/{{sample}}.{{level}}.filt.vcf.gz.tbi"),
    params:
        drop=lambda w: (
            "bcftools annotate -x " + ",".join(f"FILTER/{x}" for x in config["variant_filter"]["drop_caller_filters"]) + " -Ou -"
            if config["variant_filter"]["drop_caller_filters"] else "cat"
        ),
        min_alt=lambda w: LEVELS[w.level]["min_alt"],
        vf=config["variant_filter"],
    log: f"{OUT}/logs/filter_vcf/{{sample}}.{{level}}.log"
    conda: "envs/callers.yaml"
    shell:
        """
        (bcftools reheader --fai {REF}.fai {input} \
          | bcftools norm -a -m-any -f {REF} -c w - -Ou \
          | {params.drop} \
          | bcftools filter -m+ -s LowSupport \
              -e 'FMT/AD[0:1] < {params.min_alt} || FMT/DP[0] < {params.vf[min_dp]} || FMT/AF[0:0] < {params.vf[min_af]}' \
              -Oz -o {output.vcf}
        bcftools index -t {output.vcf}) 2> {log}
        """


# Background polishing (iDES): test each call against the per-site error rate seen
# in panel-of-normals consensus BAMs at the same level. No normals -> pass-through.

rule background_polish:
    input:
        vcf=f"{OUT}/vcf/{{sample}}.{{level}}.filt.vcf.gz",
        tbi=f"{OUT}/vcf/{{sample}}.{{level}}.filt.vcf.gz.tbi",
        normals=expand(f"{OUT}/consensus/{{n}}.{{{{level}}}}.bam", n=NORMALS),
        fai=ancient(REF + ".fai"),
    output:
        vcf=f"{OUT}/vcf/{{sample}}.{{level}}.vcf.gz",
        tbi=f"{OUT}/vcf/{{sample}}.{{level}}.vcf.gz.tbi",
    params:
        bg=config["background"],
        pile=f"{OUT}/vcf/{{sample}}.{{level}}.pon_pileup.vcf.gz",
    threads: 2
    resources: mem_mb=4000
    log: f"{OUT}/logs/background_polish/{{sample}}.{{level}}.log"
    conda: "envs/py.yaml"
    shell:
        """
        (if [ -z "{input.normals}" ]; then
            cp {input.vcf} {output.vcf}; cp {input.tbi} {output.tbi}
        else
            bcftools query -f '%CHROM\t%POS\n' {input.vcf} > {params.pile}.sites
            bcftools mpileup -f {REF} -R {params.pile}.sites -a FORMAT/AD,FORMAT/DP \
              -d 1000000 -q {params.bg[min_mq]} -Q {params.bg[min_bq]} -B --threads {threads} \
              {input.normals} -Ou \
            | bcftools norm -m-any -f {REF} -Oz -o {params.pile}
            python scripts/background_filter.py --calls {input.vcf} --pileup {params.pile} \
              --alpha {params.bg[alpha]} --min-normal-depth {params.bg[min_normal_depth]} \
              --error-floor {params.bg[error_floor]} --max-normal-af {params.bg[max_normal_af]} \
            | bcftools view -Oz -o {output.vcf}
            bcftools index -t {output.vcf}
            rm -f {params.pile}.sites
        fi) 2> {log}
        """


# Tumour-informed mode (Signatera / INVAR): genotype the truth loci directly on the
# consensus BAM, then aggregate across loci. Independent of the de novo caller.

rule genotype_truth:
    input:
        bam=f"{OUT}/consensus/{{sample}}.{{level}}.bam",
        truth=lambda w: f"{OUT}/truth/{samples.loc[w.sample, 'truth']}.norm.vcf.gz",
        bed=f"{OUT}/ref/panel.bed",
        normals=expand(f"{OUT}/consensus/{{n}}.{{{{level}}}}.bam", n=NORMALS),
        excl=lambda w: f"{OUT}/truth/{samples.loc[w.sample, 'truth']}.in_normals.tsv" if NORMALS else [],
    output:
        pile=temp(f"{OUT}/informed/{{sample}}.{{level}}.pileup.vcf.gz"),
        loci=f"{OUT}/informed/{{sample}}.{{level}}.loci.tsv",
        summary=f"{OUT}/informed/{{sample}}.{{level}}.summary.tsv",
    params:
        bg=config["background"],
        inf=config["informed"],
        scale=lambda w: samples.loc[w.sample, "vaf_scale"] or "1",
        ts=lambda w: config["truth_sets"][samples.loc[w.sample, "truth"]],
        sites=f"{OUT}/informed/{{sample}}.{{level}}.sites",
    threads: 2
    resources: mem_mb=4000
    log: f"{OUT}/logs/genotype_truth/{{sample}}.{{level}}.log"
    conda: "envs/py.yaml"
    shell:
        """
        (bcftools view -R {input.bed} {input.truth} | bcftools query -f '%CHROM\t%POS\n' > {params.sites}
        bcftools mpileup -f {REF} -R {params.sites} -a FORMAT/AD,FORMAT/DP -d 1000000 \
          -q {params.bg[min_mq]} -Q {params.bg[min_bq]} -B --threads {threads} \
          {input.bam} {input.normals} -Ou \
        | bcftools norm -a -m-any -f {REF} -Oz -o {output.pile}
        python scripts/informed_detect.py --pileup {output.pile} --truth {input.truth} \
          --sample {wildcards.sample} --level {wildcards.level} --vaf-scale {params.scale} \
          --truth-af-tag '{params.ts[af_tag]}' --default-truth-vaf '{params.ts[default_vaf]}' \
          --min-alt {params.inf[min_alt_molecules]} --min-loci {params.inf[min_loci_detected]} \
          --error-floor {params.bg[error_floor]} --alpha {params.inf[alpha]} \
          --exclude '{input.excl}' --bed {input.bed} --locus-alpha {params.inf[alpha]} \
          --loci-out {output.loci} --summary-out {output.summary}
        rm -f {params.sites}) 2> {log}
        """


# Benchmarking

rule truth_in_normals:
    # truth alleles also carried by the background DNA (Sample B germline): excluded
    # from scoring, since in the dilutions they are not at VAF(A) x fraction.
    # Uses the ss1 normals (deepest level); the list is the same for every level.
    input:
        truth=f"{OUT}/truth/{{truth}}.norm.vcf.gz",
        tbi=f"{OUT}/truth/{{truth}}.norm.vcf.gz.tbi",
        bed=f"{OUT}/ref/panel.bed",
        normals=expand(f"{OUT}/consensus/{{n}}.ss1.bam", n=NORMALS),
    output: f"{OUT}/truth/{{truth}}.in_normals.tsv"
    params:
        bg=config["background"],
        tmp=f"{OUT}/truth/{{truth}}.in_normals",
    threads: 2
    log: f"{OUT}/logs/truth_in_normals/{{truth}}.log"
    conda: "envs/py.yaml"
    shell:
        """
        (bcftools view -R {input.bed} {input.truth} | bcftools query -f '%CHROM\t%POS\n' > {params.tmp}.sites
        bcftools mpileup -f {REF} -R {params.tmp}.sites -a FORMAT/AD,FORMAT/DP -d 1000000 \
          -q {params.bg[min_mq]} -Q {params.bg[min_bq]} -B --threads {threads} {input.normals} -Ou \
        | bcftools norm -a -m-any -f {REF} -Oz -o {params.tmp}.pileup.vcf.gz
        python scripts/truth_in_normals.py --pileup {params.tmp}.pileup.vcf.gz --truth {input.truth} \
          --max-normal-af {params.bg[max_normal_af]} --min-normal-depth {params.bg[min_normal_depth]} \
          --out {output}
        rm -f {params.tmp}.sites {params.tmp}.pileup.vcf.gz) 2> {log}
        """


rule norm_truth:
    input: lambda w: config["truth_sets"][w.truth]["vcf"]
    # indexed: genotype_truth/benchmark query it by region (bcftools view -R)
    output:
        vcf=f"{OUT}/truth/{{truth}}.norm.vcf.gz",
        tbi=f"{OUT}/truth/{{truth}}.norm.vcf.gz.tbi",
    log: f"{OUT}/logs/norm_truth/{{truth}}.log"
    conda: "envs/callers.yaml"
    shell:
        "(bcftools norm -a -m-any -f {REF} -c w {input} -Oz -o {output.vcf} "
        "&& bcftools index -t {output.vcf}) 2> {log}"


rule benchmark:
    input:
        calls=f"{OUT}/vcf/{{sample}}.{{level}}.vcf.gz",
        truth=lambda w: f"{OUT}/truth/{samples.loc[w.sample, 'truth']}.norm.vcf.gz",
        bed=f"{OUT}/ref/panel.bed",
        excl=lambda w: f"{OUT}/truth/{samples.loc[w.sample, 'truth']}.in_normals.tsv" if NORMALS else [],
    output: f"{OUT}/benchmark/{{sample}}.{{level}}.tsv"
    params:
        ts=lambda w: config["truth_sets"][samples.loc[w.sample, "truth"]],
        scale=lambda w: samples.loc[w.sample, "vaf_scale"] or "1",
        neg=f"--negative-bed {config['negative_bed']}" if config.get("negative_bed") else "",
        bins=",".join(str(b) for b in config["benchmark"]["vaf_bins"]),
    log: f"{OUT}/logs/benchmark/{{sample}}.{{level}}.log"
    conda: "envs/py.yaml"
    shell:
        """
        python scripts/benchmark.py --calls {input.calls} --truth {input.truth} \
          --bed {input.bed} {params.neg} --vaf-scale {params.scale} \
          --truth-af-tag '{params.ts[af_tag]}' --default-truth-vaf '{params.ts[default_vaf]}' \
          --bins {params.bins} --sample {wildcards.sample} --level {wildcards.level} \
          --exclude '{input.excl}' --out {output} 2> {log}
        """


rule summarise:
    input:
        bench=expand(f"{OUT}/benchmark/{{sample}}.{{level}}.tsv", sample=BENCH_SAMPLES, level=LEVELS),
        fam=expand(f"{OUT}/metrics/{{sample}}.family_sizes.txt", sample=SAMPLES) if MODE == "duplex" else [],
        yld=expand(f"{OUT}/metrics/{{sample}}.duplex_yield_metrics.txt", sample=SAMPLES) if MODE == "duplex" else [],
        informed=expand(f"{OUT}/informed/{{sample}}.{{level}}.summary.tsv", sample=BENCH_SAMPLES, level=LEVELS),
        loci=expand(f"{OUT}/informed/{{sample}}.{{level}}.loci.tsv", sample=BENCH_SAMPLES, level=LEVELS),
    params:
        bins=",".join(str(b) for b in config["benchmark"]["vaf_bins"]),
    output:
        f"{OUT}/benchmark/summary.tsv",
        directory(f"{OUT}/benchmark/plots"),
    log: f"{OUT}/logs/summarise.log"
    conda: "envs/py.yaml"
    shell:
        """
        python scripts/summarise.py --bench {input.bench} --family-sizes {input.fam} \
          --duplex-yield {input.yld} --informed {input.informed} --informed-loci {input.loci} \
          --vaf-bins {params.bins} --outdir {OUT}/benchmark 2> {log}
        """


rule multiqc:
    input:
        expand(f"{OUT}/qc/fastqc/{{sample}}", sample=SAMPLES) if config.get("run_fastqc", True) else [],
        # resume_from_bams: use only the fastp reports that exist, so a missing report
        # never pulls the download/alignment chain back in
        [f for f in expand(f"{OUT}/qc/fastp/{{sample}}.fastp.json", sample=SAMPLES)
         if not config.get("resume_from_bams", False) or os.path.exists(f)],
        expand(f"{OUT}/metrics/{{sample}}.family_size_histogram.txt", sample=SAMPLES),
        expand(f"{OUT}/qc/mosdepth/{{sample}}.{{level}}.mosdepth.summary.txt", sample=SAMPLES, level=LEVELS),
        expand(f"{OUT}/qc/fragments/{{sample}}.{{level}}.stats.txt", sample=SAMPLES, level=LEVELS),
    output: f"{OUT}/qc/multiqc_report.html"
    log: f"{OUT}/logs/multiqc.log"
    conda: "envs/core.yaml"
    shell: "multiqc -f -o {OUT}/qc {OUT}/qc {OUT}/metrics &> {log}"
