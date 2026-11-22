#!/usr/bin/env python3
"""Tumour-informed detection at known truth loci (Signatera / INVAR style).

Pileup column 0 is the test sample; any further columns are panel-of-normals and give
locus-specific background rates (otherwise the error floor is used).

Per locus:    alt molecules, depth, VAF, detected = alt >= min_alt AND the ALT count is
              significantly above that locus's panel-of-normals background
              (binomial p <= --locus-alpha), so single-strand error noise at deep loci
              does not count as a detected locus
Per sample:   fraction of loci detected, Signatera call (>= min_loci detected),
              INVAR-style aggregate call (Poisson test of total ALT across all loci
              against sum(depth_i * background_i)), mean observed vs expected VAF.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from vcfutil import binom_sf, in_bed_set, poisson_sf, read_bed_set, read_excluded, read_pileup_ad, read_truth


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pileup", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--vaf-scale", type=float, default=1.0)
    ap.add_argument("--truth-af-tag", default="")
    ap.add_argument("--default-truth-vaf", default="")
    ap.add_argument("--min-alt", type=int, default=2)
    ap.add_argument("--min-loci", type=int, default=2)
    ap.add_argument("--error-floor", type=float, default=5e-5)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--loci-out", required=True)
    ap.add_argument("--summary-out", required=True)
    ap.add_argument("--exclude", default="", help="truth_in_normals.tsv: loci shared with the normals")
    ap.add_argument("--bed", default="", help="panel BED: only truth loci inside it are scored")
    ap.add_argument("--locus-alpha", type=float, default=0.01)
    a = ap.parse_args()

    dv = float(a.default_truth_vaf) if a.default_truth_vaf else None
    truth = read_truth(a.truth, a.truth_af_tag, dv)
    excl = read_excluded(a.exclude)
    bed = read_bed_set(a.bed) if a.bed else None
    truth = {k: v for k, v in truth.items()
             if k not in excl and (bed is None or in_bed_set(bed, k[0], k[1]))}
    pile = read_pileup_ad(a.pileup)

    n_loci = n_det = n_cov = 0
    tot_alt = 0
    lam = 0.0
    vafs, exp_vafs = [], []
    with open(a.loci_out, "w") as o:
        o.write("sample\tlevel\tchrom\tpos\tref\talt\texpected_vaf\tdepth\talt_molecules\tvaf\tbackground_rate\tlocus_p\tdetected\n")
        for key, tv in truth.items():
            n_loci += 1
            per = pile.get(key, [(0, 0)])
            alt, dp = per[0]
            bg_alt = sum(x[0] for x in per[1:])
            bg_dp = sum(x[1] for x in per[1:])
            bg = max(a.error_floor, (bg_alt + 0.5) / (bg_dp + 1.0)) if bg_dp >= 50 else a.error_floor
            exp_v = tv * a.vaf_scale if tv is not None else None
            p_loc = binom_sf(alt, dp, bg) if dp else 1.0
            det = int(alt >= a.min_alt and p_loc <= a.locus_alpha)
            if dp > 0:
                n_cov += 1
                vafs.append(alt / dp)
                if exp_v is not None:
                    exp_vafs.append(exp_v)
            n_det += det
            tot_alt += alt
            lam += dp * bg
            o.write("\t".join(str(x) for x in (
                a.sample, a.level, key[0], key[1], key[2], key[3],
                "" if exp_v is None else f"{exp_v:.5f}", dp, alt,
                f"{alt / dp:.6f}" if dp else "", f"{bg:.2e}", f"{p_loc:.3g}", det)) + "\n")

    p_agg = poisson_sf(tot_alt, lam) if n_cov else 1.0
    with open(a.summary_out, "w") as o:
        o.write("sample\tlevel\tn_truth_loci\tn_covered\tn_detected\tfrac_detected\ttotal_alt_molecules\t"
                "expected_background_alt\taggregate_p\tmean_observed_vaf\tmean_expected_vaf\t"
                "call_signatera\tcall_invar\n")
        o.write("\t".join(str(x) for x in (
            a.sample, a.level, n_loci, n_cov, n_det,
            f"{n_det / n_cov:.4f}" if n_cov else "", tot_alt, f"{lam:.2f}", f"{p_agg:.3g}",
            f"{sum(vafs) / len(vafs):.6f}" if vafs else "",
            f"{sum(exp_vafs) / len(exp_vafs):.6f}" if exp_vafs else "",
            int(n_det >= a.min_loci), int(p_agg <= a.alpha))) + "\n")


if __name__ == "__main__":
    main()
