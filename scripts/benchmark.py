#!/usr/bin/env python3
"""Benchmark normalised calls against a normalised truth VCF within the panel.

Sensitivity is binned by EXPECTED VAF (truth VAF x dilution scale), because
missed variants have no observed VAF. PPV is binned by OBSERVED VAF, because
false positives have no expected VAF. If a known-negative BED is given, only
non-truth calls inside it count as FP; others are reported as 'unclassified'.
"""
import argparse
import bisect
import gzip
import sys
from collections import defaultdict


def open_text(path):
    # gzip detected by magic bytes, so .bed.gz, .bgz or a misnamed file all work
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if gz else open(path)


def read_bed(path):
    iv = defaultdict(list)
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            iv[c].append((int(s), int(e)))
    merged, size = {}, 0
    for c, lst in iv.items():
        lst.sort()
        out = []
        for s, e in lst:
            if out and s <= out[-1][1]:
                out[-1][1] = max(out[-1][1], e)
            else:
                out.append([s, e])
        merged[c] = ([x[0] for x in out], [x[1] for x in out])
        size += sum(e - s for s, e in out)
    return merged, size


def overlap_size(a, b):
    """Total bases shared by two merged interval sets from read_bed."""
    tot = 0
    for c, (sa, ea) in a.items():
        if c not in b:
            continue
        sb, eb = b[c]
        i = j = 0
        while i < len(sa) and j < len(sb):
            lo, hi = max(sa[i], sb[j]), min(ea[i], eb[j])
            if lo < hi:
                tot += hi - lo
            if ea[i] < eb[j]:
                i += 1
            else:
                j += 1
    return tot


def intersect(a, b):
    """Intersection of two merged interval sets from read_bed, same structure."""
    out = {}
    for c, (sa, ea) in a.items():
        if c not in b:
            continue
        sb, eb = b[c]
        i = j = 0
        st, en = [], []
        while i < len(sa) and j < len(sb):
            lo, hi = max(sa[i], sb[j]), min(ea[i], eb[j])
            if lo < hi:
                st.append(lo); en.append(hi)
            if ea[i] < eb[j]:
                i += 1
            else:
                j += 1
        if st:
            out[c] = (st, en)
    return out


def masked_negative_bp(neg, bed, mask):
    """Known-negative panel bases that fall inside the mask."""
    return overlap_size(intersect(neg, bed), mask)


def in_bed(bed, chrom, pos1):
    if chrom not in bed:
        return False
    starts, ends = bed[chrom]
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    return i >= 0 and pos1 - 1 < ends[i]


def info_get(info, tag):
    for kv in info.split(";"):
        if kv.startswith(tag + "="):
            return kv.split("=", 1)[1].split(",")[0]
    return None


def read_vcf(path, pass_only, af_tag=None):
    op = gzip.open if path.endswith(".gz") else open
    out = {}
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, flt, info = f[0], int(f[1]), f[3], f[4], f[6], f[7]
            if alt in (".", "*") or alt.startswith("<"):
                continue
            if pass_only and flt not in ("PASS", "."):
                continue
            af = None
            if af_tag:
                v = info_get(info, af_tag)
                af = float(v) if v not in (None, ".") else None
            elif len(f) > 9:
                fmt = f[8].split(":")
                if "AF" in fmt:
                    v = f[9].split(":")[fmt.index("AF")].split(",")[0]
                    af = float(v) if v != "." else None
                if af is None:
                    v = info_get(info, "AF")
                    af = float(v) if v not in (None, ".") else None
            out[(chrom, pos, ref, alt)] = af
    return out


def bin_label(v, bins):
    if v is None:
        return "NA"
    for lo, hi in zip(bins[:-1], bins[1:]):
        if lo <= v < hi or (hi == bins[-1] and v == hi):
            return f"[{lo:g},{hi:g})"
    return "NA"


COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def sub_class(ref, alt):
    """Pyrimidine-centred substitution class (C>T covers G>A etc.); 'indel' otherwise."""
    if len(ref) != 1 or len(alt) != 1:
        return "indel"
    if ref in "AG":
        ref, alt = COMP[ref], COMP[alt]
    return f"{ref}>{alt}"


def safe(n, d):
    return round(n / d, 4) if d else ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bed", required=True)
    ap.add_argument("--negative-bed")
    ap.add_argument("--vaf-scale", type=float, default=1.0)
    ap.add_argument("--truth-af-tag", default="")
    ap.add_argument("--default-truth-vaf", default="")
    ap.add_argument("--bins", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--exclude", default="", help="truth_in_normals.tsv: truth alleles shared with the normals")
    ap.add_argument("--mask-bed", default="", help="BED of sites to ignore entirely (e.g. spike-in standards); "
                    "calls there are dropped and the bases leave the FP space")
    a = ap.parse_args()

    bins = [float(x) for x in a.bins.split(",")]
    bed, bed_size = read_bed(a.bed)
    neg, neg_size = read_bed(a.negative_bed) if a.negative_bed else (None, 0)
    default_vaf = float(a.default_truth_vaf) if a.default_truth_vaf else None

    truth_raw = read_vcf(a.truth, pass_only=False, af_tag=a.truth_af_tag or "__none__")
    excl = set()
    if a.exclude:
        with open(a.exclude) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                excl.add((f[0], int(f[1]), f[2], f[3]))
    truth = {}
    n_shared = 0
    for k, af in truth_raw.items():
        if not in_bed(bed, k[0], k[1]):
            continue
        if k in excl:
            n_shared += 1
            continue
        base = af if af is not None else default_vaf
        truth[k] = base * a.vaf_scale if base is not None else None

    mask, _ = read_bed(a.mask_bed) if a.mask_bed else (None, 0)
    calls = {k: v for k, v in read_vcf(a.calls, pass_only=True).items()
             if in_bed(bed, k[0], k[1]) and not (mask is not None and in_bed(mask, k[0], k[1]))}
    n_masked_truth = 0
    if mask is not None:
        for k in [k for k in truth if in_bed(mask, k[0], k[1])]:
            del truth[k]; n_masked_truth += 1

    tp = set(calls) & set(truth)
    fn = set(truth) - set(calls)
    other = set(calls) - set(truth)
    if neg is not None:
        fp = {k for k in other if in_bed(neg, k[0], k[1])}
        unclassified = other - fp
        # FP rate denominator: known-negative bases INSIDE the panel (calls only come from
        # the panel; the genome-wide negative BED is ~10 Mb and would deflate the rate ~100x)
        fp_space = overlap_size(neg, bed)
        if mask is not None:
            fp_space -= masked_negative_bp(neg, bed, mask)
    else:
        fp, unclassified, fp_space = other, set(), bed_size

    rows = []
    sens = defaultdict(lambda: [0, 0])
    for k in tp:
        sens[bin_label(truth[k], bins)][0] += 1
    for k in fn:
        sens[bin_label(truth[k], bins)][1] += 1
    ppv = defaultdict(lambda: [0, 0])
    for k in tp:
        ppv[bin_label(calls[k], bins)][0] += 1
    for k in fp:
        ppv[bin_label(calls[k], bins)][1] += 1

    order = [bin_label((lo + hi) / 2, bins) for lo, hi in zip(bins[:-1], bins[1:])] + ["NA"]
    # vaf_scale 0 = negative control (Sample B): no expected positives, report FP/PPV only
    if a.vaf_scale == 0:
        sens.clear()
        fn = set()
    for b in order:
        if b in sens:
            t, f = sens[b]
            rows.append(("sensitivity", b, t, f, "", safe(t, t + f)))
    for b in order:
        if b in ppv:
            t, f = ppv[b]
            rows.append(("ppv", b, t, "", f, safe(t, t + f)))
    if a.vaf_scale != 0:
        rows.append(("sensitivity", "all", len(tp), len(fn), "", safe(len(tp), len(tp) + len(fn))))
    rows.append(("ppv", "all", len(tp), "", len(fp), safe(len(tp), len(tp) + len(fp))))
    rows.append(("fp_per_mb", "negative" if neg is not None else "target", "", "", len(fp),
                 round(len(fp) / (fp_space / 1e6), 3) if fp_space else ""))
    # residual error spectrum (iDES: C>T deamination / G>T oxidation dominate if polishing is missing)
    spec = defaultdict(int)
    for k in fp:
        spec[sub_class(k[2], k[3])] += 1
    for cls in ("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "indel"):
        rows.append(("fp_spectrum", cls, "", "", spec[cls], safe(spec[cls], len(fp))))
    rows.append(("unclassified_calls", "target", "", "", "", len(unclassified)))
    rows.append(("fp_space_bp", "negative_in_panel" if neg is not None else "target", "", "", "", fp_space))
    rows.append(("truth_in_target", "all", "", "", "", len(truth)))
    rows.append(("truth_shared_with_normals", "excluded", "", "", "", n_shared))
    if mask is not None:
        rows.append(("masked_truth", "excluded", "", "", "", n_masked_truth))

    with open(a.out, "w") as o:
        o.write("sample\tlevel\tmetric\tbin\ttp\tfn\tfp\tvalue\n")
        for r in rows:
            o.write("\t".join(str(x) for x in (a.sample, a.level) + r) + "\n")
    print(f"{a.sample} {a.level}: TP={len(tp)} FN={len(fn)} FP={len(fp)}", file=sys.stderr)


if __name__ == "__main__":
    main()
