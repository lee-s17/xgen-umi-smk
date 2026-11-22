#!/usr/bin/env python3
"""iDES-style background polishing.

For each call, pool ALT/depth across the panel-of-normals pileup at that site and
allele, estimate the site error rate (with pseudocount and a floor), and test whether
the sample's ALT count is significantly above it (binomial). Adds INFO BG_AF, BG_DP,
BG_P and FILTER Background. Writes VCF to stdout.

A call whose allele is plainly present in the normals (pooled AF >= --max-normal-af) is
flagged PoN regardless of the test: that is germline variation of the background DNA or a
synthetic spike-in shared by every library (e.g. SEQC2 Accugenomics internal standards),
not a somatic variant. A binomial test alone keeps these when the sample is deep.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vcfutil import binom_sf, read_pileup_ad, vopen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--pileup", required=True)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--min-normal-depth", type=int, default=50)
    ap.add_argument("--error-floor", type=float, default=5e-5)
    ap.add_argument("--max-normal-af", type=float, default=0.01)
    a = ap.parse_args()

    pile = read_pileup_ad(a.pileup)
    n_flag = n_tot = n_pon = 0
    with vopen(a.calls) as fh:
        for line in fh:
            if line.startswith("##"):
                sys.stdout.write(line)
                continue
            if line.startswith("#CHROM"):
                sys.stdout.write('##INFO=<ID=BG_AF,Number=1,Type=Float,Description="Pooled panel-of-normals ALT fraction at this site/allele">\n')
                sys.stdout.write('##INFO=<ID=BG_DP,Number=1,Type=Integer,Description="Pooled panel-of-normals depth at this site">\n')
                sys.stdout.write('##INFO=<ID=BG_P,Number=1,Type=Float,Description="Binomial P(ALT >= observed | background rate)">\n')
                sys.stdout.write(f'##FILTER=<ID=Background,Description="Not above panel-of-normals background (p > {a.alpha})">\n')
                sys.stdout.write(f'##FILTER=<ID=PoN,Description="Allele present in panel of normals at AF >= {a.max_normal_af} (germline or spike-in)">\n')
                sys.stdout.write(line)
                continue
            f = line.rstrip("\n").split("\t")
            key = (f[0], int(f[1]), f[3], f[4])
            fmt = f[8].split(":")
            vals = f[9].split(":")
            ad = vals[fmt.index("AD")].split(",") if "AD" in fmt else ["0", "0"]
            alt = int(ad[1]) if len(ad) > 1 and ad[1] != "." else 0
            dp = int(vals[fmt.index("DP")]) if "DP" in fmt and vals[fmt.index("DP")] != "." else sum(int(x) for x in ad if x != ".")

            bg_alt = sum(x[0] for x in pile.get(key, []))
            bg_dp = sum(x[1] for x in pile.get(key, []))
            if bg_dp >= a.min_normal_depth:
                rate = max(a.error_floor, (bg_alt + 0.5) / (bg_dp + 1.0))
            else:
                rate = a.error_floor
            p = binom_sf(alt, dp, rate) if dp else 1.0
            bg_af = bg_alt / bg_dp if bg_dp else 0.0
            f[7] = (f[7] + ";" if f[7] != "." else "") + f"BG_AF={bg_af:.6f};BG_DP={bg_dp};BG_P={p:.3g}"
            n_tot += 1
            flags = []
            if bg_dp >= a.min_normal_depth and bg_af >= a.max_normal_af:
                flags.append("PoN")
                n_pon += 1
            if p > a.alpha:
                flags.append("Background")
            if flags:
                n_flag += 1
                f[6] = ";".join(flags) if f[6] in ("PASS", ".") else f[6] + ";" + ";".join(flags)
            sys.stdout.write("\t".join(f) + "\n")
    print(f"background_filter: flagged {n_flag}/{n_tot} calls ({n_pon} present in normals)", file=sys.stderr)


if __name__ == "__main__":
    main()
