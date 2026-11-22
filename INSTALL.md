# Installation

The workflow needs one small environment to run Snakemake. Snakemake then builds the
three tool environments in `envs/` itself, so the tool versions are pinned by the repo.

## 1. Snakemake environment

```bash
# once per machine: bioconda's recommended channel set-up
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

mamba create -n xgen-umi -c conda-forge -c bioconda "snakemake>=9" "python>=3.11" pandas
mamba activate xgen-umi
snakemake --version
```

## 2. Tool environments (built by Snakemake)

```bash
git clone https://github.com/lee-s17/xgen-umi-smk && cd xgen-umi-smk
snakemake --profile profiles/local --conda-create-envs-only   # builds envs/*.yaml once
```

| Env file | Tools |
| --- | --- |
| `envs/core.yaml` | fgbio 4.x (Java 17), bwa-mem2, samtools, fastp, FastQC, MultiQC, mosdepth, bedtools, sra-tools, pigz |
| `envs/callers.yaml` | VarDict-java, R, GATK4 (Mutect2), bcftools |
| `envs/py.yaml` | Python, pandas, matplotlib, bcftools |

Sanity checks inside the core env (`conda activate` the path Snakemake prints under `.snakemake/conda/`):

```bash
fgbio --version 2>&1 | head -1                    # 4.x; needs Java 17
java -version 2>&1 | head -1                      # 17
bwa-mem2 version                                  # the binary may print 2.2.1 even for the 2.3 package
samtools sort --help 2>&1 | grep template-coordinate
fastqc --version
```

In the callers env, the VarDict command is `vardict-java` (not `vardict`).

## 3. bwa-mem2 index (one-off)

```bash
REF=resources/ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna   # or set ref: in config.yaml
samtools faidx "$REF"
bwa-mem2 index "$REF"            # ~1 h, needs ~60-70 GB RAM, writes ~17 GB next to the FASTA
ls "$REF".*                      # .0123 .amb .ann .bwt.2bit.64 .pac
```

Or `snakemake --profile profiles/local bwa_index_all`. The alignment rules stop with a clear
message if the index is missing, so a full run never starts a rebuild by accident. Build it on
a machine with enough RAM and copy the five index files over if needed.

## 4. SRA toolkit first run

```bash
vdb-config --interactive       # set a download cache on a disk with space, then exit
prefetch --version
```

## 5. Reference, panel and truth checks

```bash
head -1 "$REF" | cut -c1-30                                    # >chr1 (no-alt set uses chr names)
head -3 resources/bed/LBx_IDT_liftover_hg38.bed               # must also be chr-prefixed
bcftools view -h resources/truth/KnownPositives_hg19ToHg38.vcf.gz | grep -m1 contig
```

The SEQC2 truth files and IDT panel BED come from figshare (see README). Lift hg19 files to
hg38 first (`mamba install -c bioconda ucsc-liftover` or `crossmap`, chain hg19ToHg38).

## 6. Run

```bash
snakemake -n --profile profiles/local     # dry run
snakemake --profile profiles/local
```

`profiles/local/config.yaml` sets cores, a total memory cap (`mem_mb`; set to total RAM minus
~12 GB) and conda deployment. Run long jobs inside `tmux` or `screen`.
