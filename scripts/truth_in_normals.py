#!/usr/bin/env python3
"""List truth alleles that are present in the panel of normals.

SEQC2 Sample A is a pool of cell lines; some of its "known positives" are common
germline variants that the background DNA (Sample B) also carries. In a dilution of
A into B those sit at B's genotype VAF, not at VAF(A) x fraction, and the PoN filter
correctly removes them as germline. Scoring them as missed somatic variants would
wrongly lower sensitivity in the high-VAF bins, so they are listed here and excluded
from the benchmark and from tumour-informed detection.

Input: normalised bcftools mpileup of the normal BAMs at truth sites.
Output TSV: chrom, pos, ref, alt, normal_alt, normal_depth, normal_af
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vcfutil import read_pileup_ad, read_truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pileup", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--max-normal-af", type=float, default=0.01)
    ap.add_argument("--min-normal-depth", type=int, default=50)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    truth = read_truth(a.truth)
    pile = read_pileup_ad(a.pileup)
    n = 0
    with open(a.out, "w") as o:
        o.write("chrom\tpos\tref\talt\tnormal_alt\tnormal_depth\tnormal_af\n")
        for key, per in pile.items():
            if key not in truth:
                continue
            alt = sum(x[0] for x in per)
            dp = sum(x[1] for x in per)
            if dp >= a.min_normal_depth and alt / dp >= a.max_normal_af:
                n += 1
                o.write("\t".join(str(x) for x in (*key, alt, dp, f"{alt / dp:.4f}")) + "\n")
    print(f"truth alleles present in normals: {n}", file=sys.stderr)


if __name__ == "__main__":
    main()
