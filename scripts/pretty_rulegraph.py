#!/usr/bin/env python3
"""Restyle `snakemake --rulegraph` output: stage clusters, tool sublabels,
left-to-right layout, reference edges dimmed, the `all` target dropped.

usage: snakemake --rulegraph | python scripts/pretty_rulegraph.py [--show-ref] | dot -Tsvg > rulegraph.svg

Reference edges (dict, fai, BED) feed nearly every step and tangle the graph,
so by default they are hidden and the reference box is pinned beside Input.
"""
import re
import sys

# stage -> (title, fill, border, rules)
STAGES = [
    ("input", "Input", "#EEF4FB", "#3B6EA8",
     ["fetch_sra", "fastp", "fastqc"]),
    ("ref", "Reference and targets", "#F3F3F1", "#8A8A85",
     ["bwa_index", "bwa_index_all", "ref_dict", "ref_fai", "prep_bed"]),
    ("umi", "UMI consensus (fgbio)", "#EFEEFD", "#4B3FA8",
     ["fastq_to_ubam", "align_raw", "group_umis", "call_consensus",
      "align_consensus", "filter_clip"]),
    ("call", "Variant calling and filtering", "#FDF1EA", "#B0562A",
     ["vardict", "filter_vcf", "background_polish"]),
    ("qc", "QC", "#EDF7EF", "#2F7D45",
     ["duplex_metrics", "fragment_stats", "mosdepth", "multiqc"]),
    ("bench", "Truth and benchmarking", "#FBF4E3", "#9A6B12",
     ["norm_truth", "truth_in_normals", "genotype_truth", "benchmark",
      "summarise"]),
]

TOOLS = {
    "fetch_sra": "fasterq-dump",
    "fastp": "adapter/quality trim",
    "fastqc": "FastQC",
    "ref_dict": "samtools dict",
    "ref_fai": "samtools faidx",
    "bwa_index": "bwa-mem2 index",
    "prep_bed": "panel BED",
    "fastq_to_ubam": "FastqToBam 3M2S146T",
    "align_raw": "bwa-mem2 + ZipperBams",
    "group_umis": "GroupReadsByUmi paired",
    "call_consensus": "CallDuplexConsensusReads",
    "align_consensus": "bwa-mem2 + ZipperBams",
    "filter_clip": "FilterConsensusReads, ClipBam\nx4 levels: ss1 ss3 duplex strict",
    "vardict": "VarDict-java",
    "filter_vcf": "bcftools norm/filter",
    "background_polish": "PoN + binomial background",
    "duplex_metrics": "CollectDuplexSeqMetrics",
    "fragment_stats": "insert sizes",
    "mosdepth": "consensus depth",
    "multiqc": "MultiQC report",
    "norm_truth": "bcftools norm",
    "truth_in_normals": "truth alleles in Sample B",
    "genotype_truth": "tumour-informed test",
    "benchmark": "sens / PPV / FP per Mb",
    "summarise": "tables + plots",
}

REF_RULES = {"ref_dict", "ref_fai", "bwa_index", "prep_bed"}


def parse(text):
    labels = dict(re.findall(r'^\s*(\d+)\s*\[.*?label\s*=\s*"([^"]+)"', text, re.M))
    edges = re.findall(r'^\s*(\d+)\s*->\s*(\d+)', text, re.M)
    return labels, [(labels[a], labels[b]) for a, b in edges]


def main():
    show_ref = "--show-ref" in sys.argv[1:]
    labels, edges = parse(sys.stdin.read())
    rules = set(labels.values()) - {"all"}
    edges = [(a, b) for a, b in edges if "all" not in (a, b)]
    placed = set()

    out = ['digraph rulegraph {',
           '  graph [rankdir=LR, bgcolor=white, nodesep=0.25, ranksep=0.45,',
           '         fontname="Helvetica", fontsize=11, newrank=true];',
           '  node  [shape=box, style="rounded,filled", fontname="Helvetica",',
           '         fontsize=10, penwidth=1.4, margin="0.12,0.06"];',
           '  edge  [color="#7A7A75", arrowsize=0.6, penwidth=1.1];']
    for key, title, fill, border, members in STAGES:
        present = [r for r in members if r in rules]
        if not present:
            continue
        out.append(f'  subgraph cluster_{key} {{')
        out.append(f'    label="{title}"; labeljust=l; fontcolor="{border}";')
        out.append(f'    style="rounded"; color="{border}"; penwidth=1.0;')
        for r in present:
            tool = TOOLS.get(r, "").replace("\n", "<br/>")
            lab = f'<<b>{r}</b>' + (f'<br/><font point-size="8">{tool}</font>' if tool else "") + '>'
            out.append(f'    "{r}" [label={lab}, fillcolor="{fill}", color="{border}"];')
            placed.add(r)
        out.append('  }')
    for r in sorted(rules - placed):
        out.append(f'  "{r}";')
    for a, b in edges:
        if a in REF_RULES:
            if show_ref:
                out.append(f'  "{a}" -> "{b}" [style=dashed, color="#BDBDB8"];')
        else:
            out.append(f'  "{a}" -> "{b}";')
    if not show_ref:
        # keep the reference box in the first column, under Input
        ref = [r for r in ("ref_dict", "ref_fai", "prep_bed") if r in rules]
        if ref and "fetch_sra" in rules:
            out.append('  { rank=same; "fetch_sra"; ' + " ".join(f'"{r}";' for r in ref) + ' }')
        out.append('  note [shape=plaintext, style="", fontsize=9, fontcolor="#5F5E5A", '
                   'label="Reference edges hidden: dict/fai feed alignment and calling; '
                   'the panel BED feeds calling, QC and benchmarking."];')
    out.append('}')
    print("\n".join(out))


if __name__ == "__main__":
    main()
