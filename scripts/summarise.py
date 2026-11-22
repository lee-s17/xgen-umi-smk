#!/usr/bin/env python3
"""Pool per-sample benchmark tables and draw report figures (PNG)."""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SURFACE, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
    "lines.markersize": 7, "font.size": 10, "legend.frameon": False,
})


def colour_map(keys):
    # colour follows the entity in a fixed order, never its rank
    return {k: SERIES[i % len(SERIES)] for i, k in enumerate(keys)}


def save(fig, outdir, name):
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, name), dpi=200)
    plt.close(fig)


def by_bin(df, metric, levels, cmap, ylabel, title, outdir, name):
    d = df[(df.metric == metric) & (df.bin != "all") & (df.bin != "NA")]
    if d.empty:
        return
    bins = list(dict.fromkeys(d.bin))
    fig, ax = plt.subplots(figsize=(7, 4))
    for lv in levels:
        x = d[d.level == lv].set_index("bin").reindex(bins)
        ax.plot(range(len(bins)), x.value, marker="o", color=cmap[lv], label=lv)
    # n per bin (variants for sensitivity, calls for PPV) from the first level, under each tick
    n = d[d.level == levels[0]].set_index("bin").reindex(bins)
    tot = (n.tp + (n.fn if metric == "sensitivity" else n.fp)).fillna(0).astype(int)
    labels = [f"{b}\nn={t}" for b, t in zip(bins, tot)]
    ax.set_xticks(range(len(bins)), labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("VAF bin")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.legend(title="Family-size level", ncols=2)
    save(fig, outdir, name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bench", nargs="+", required=True)
    ap.add_argument("--family-sizes", nargs="*", default=[])
    ap.add_argument("--duplex-yield", nargs="*", default=[])
    ap.add_argument("--informed", nargs="*", default=[])
    ap.add_argument("--informed-loci", nargs="*", default=[])
    ap.add_argument("--vaf-bins", default="0,0.001,0.0025,0.005,0.01,0.02,0.05,1")
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    pdir = os.path.join(a.outdir, "plots")
    os.makedirs(pdir, exist_ok=True)

    raw = pd.concat([pd.read_csv(f, sep="\t") for f in a.bench], ignore_index=True)
    raw.to_csv(os.path.join(a.outdir, "per_sample.tsv"), sep="\t", index=False)
    levels = list(dict.fromkeys(raw.level))
    cmap = colour_map(levels)

    # pool counts across replicates, then recompute rates
    rates = raw[raw.metric.isin(["sensitivity", "ppv"])].copy()
    for c in ("tp", "fn", "fp"):
        rates[c] = pd.to_numeric(rates[c], errors="coerce").fillna(0)
    pooled = rates.groupby(["level", "metric", "bin"], sort=False)[["tp", "fn", "fp"]].sum().reset_index()
    denom = pooled.tp + pooled.fn.where(pooled.metric == "sensitivity", pooled.fp)
    pooled["value"] = (pooled.tp / denom).where(denom > 0).round(4)

    fpmb = raw[raw.metric == "fp_per_mb"].copy()
    fpmb["value"] = pd.to_numeric(fpmb.value, errors="coerce")
    fp_level = fpmb.groupby("level", sort=False).value.mean().reset_index()
    fp_level["metric"], fp_level["bin"] = "fp_per_mb_mean", fpmb.bin.iloc[0] if len(fpmb) else ""
    pd.concat([pooled, fp_level], ignore_index=True).to_csv(
        os.path.join(a.outdir, "summary.tsv"), sep="\t", index=False)

    by_bin(pooled, "sensitivity", levels, cmap, "Sensitivity",
           "Sensitivity by expected VAF (pooled replicates)", pdir, "sensitivity_by_vaf.png")
    by_bin(pooled, "ppv", levels, cmap, "PPV",
           "PPV by observed VAF (pooled replicates)", pdir, "ppv_by_vaf.png")

    spec = raw[raw.metric == "fp_spectrum"].copy()
    if len(spec):
        spec["fp"] = pd.to_numeric(spec.fp, errors="coerce").fillna(0)
        sp = spec.groupby(["level", "bin"], sort=False).fp.sum().unstack("bin").reindex(levels)
        sp = sp.div(sp.sum(axis=1).replace(0, 1), axis=0)
        sp.to_csv(os.path.join(a.outdir, "fp_spectrum.tsv"), sep="\t")
        fig, ax = plt.subplots(figsize=(6, 3.5))
        bottom = pd.Series(0.0, index=sp.index)
        for i, cls in enumerate(sp.columns):
            ax.bar(sp.index, sp[cls], bottom=bottom, color=SERIES[i % len(SERIES)], width=0.6,
                   label=cls, edgecolor=SURFACE, linewidth=1)
            bottom += sp[cls]
        ax.set_ylabel("Fraction of false positives")
        ax.set_title("Residual FP error spectrum by level", loc="left")
        ax.legend(ncols=4, fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.18))
        ax.grid(axis="x", visible=False)
        save(fig, pdir, "fp_spectrum.png")

    if a.informed:
        inf = pd.concat([pd.read_csv(f, sep="\t") for f in a.informed], ignore_index=True)
        inf.to_csv(os.path.join(a.outdir, "informed_summary.tsv"), sep="\t", index=False)
        # aggregate (INVAR-style) evidence per library: negative control must stay below the line
        inf["neglog10p"] = -pd.Series(pd.to_numeric(inf.aggregate_p, errors="coerce")).clip(lower=1e-300).map(
            lambda v: __import__("math").log10(v))
        samples_ = list(dict.fromkeys(inf["sample"]))
        fig, ax = plt.subplots(figsize=(8, 4))
        width = 0.8 / max(1, len(levels))
        for j, lv in enumerate(levels):
            d = inf[inf.level == lv].set_index("sample").reindex(samples_)
            ax.bar([i + j * width for i in range(len(samples_))], d.neglog10p.clip(upper=50),
                   width=width, color=cmap[lv], label=lv, edgecolor=SURFACE)
        ax.axhline(2, color=INK2, linestyle="--", linewidth=1)
        ax.annotate("p = 0.01", (-0.45, 2), va="bottom", ha="left", color=INK2, fontsize=8)
        short = [s_.replace("_IDT2_ST04_25ng", "") for s_ in samples_]
        ax.set_xticks([i + width * (len(levels) - 1) / 2 for i in range(len(samples_))], short,
                      rotation=20, ha="right", fontsize=8)
        ax.set_ylabel("-log10 aggregate p (capped at 50)")
        ax.set_title("Tumour-informed aggregate detection per library", loc="left")
        ax.legend(title="Family-size level", ncols=4, fontsize=8, loc="upper center",
                  bbox_to_anchor=(0.5, -0.28))
        ax.grid(axis="x", visible=False)
        save(fig, pdir, "informed_detection.png")

    if a.informed_loci:
        loci = pd.concat([pd.read_csv(f, sep="\t") for f in a.informed_loci], ignore_index=True)
        loci = loci[(loci.depth > 0) & loci.expected_vaf.notna() & (loci.expected_vaf > 0)]
        if len(loci):
            edges = [float(x) for x in a.vaf_bins.split(",")]
            lab = [f"[{lo:g},{hi:g})" for lo, hi in zip(edges[:-1], edges[1:])]
            loci["bin"] = pd.cut(loci.expected_vaf, edges, labels=lab, right=False)
            g = loci.groupby(["level", "bin"], observed=True).detected.agg(["mean", "size"]).reset_index()
            g.to_csv(os.path.join(a.outdir, "informed_by_vaf.tsv"), sep="\t", index=False)
            fig, ax = plt.subplots(figsize=(7, 4))
            bins_ = [b for b in lab if b in set(g.bin.astype(str))]
            for lv in levels:
                x = g[g.level == lv].assign(bin=lambda t: t.bin.astype(str)).set_index("bin").reindex(bins_)
                ax.plot(range(len(bins_)), x["mean"], marker="o", color=cmap[lv], label=lv)
            nn = g[g.level == levels[0]].assign(bin=lambda t: t.bin.astype(str)).set_index("bin").reindex(bins_)["size"]
            ax.set_xticks(range(len(bins_)), [f"{b}\nn={int(v) if v == v else 0}" for b, v in zip(bins_, nn)],
                          rotation=30, ha="right", fontsize=8)
            ax.set_ylim(0, 1.02)
            ax.set_xlabel("Expected VAF bin")
            ax.set_ylabel("Loci with >= min ALT molecules")
            ax.set_title("Tumour-informed per-locus detection by expected VAF", loc="left")
            ax.legend(title="Family-size level", ncols=2)
            save(fig, pdir, "informed_by_vaf.png")

    if len(fp_level):
        fig, ax = plt.subplots(figsize=(5, 3.5))
        ax.bar(fp_level.level, fp_level.value, color=[cmap[l] for l in fp_level.level], width=0.6)
        for i, v in enumerate(fp_level.value):
            ax.annotate(f"{v:.1f}", (i, v), ha="center", va="bottom", color=INK2, fontsize=9)
        ax.set_ylabel("False positives per Mb")
        ax.set_title("False-positive rate by level", loc="left")
        ax.grid(axis="x", visible=False)
        save(fig, pdir, "fp_per_mb.png")

    if a.family_sizes:
        fig, ax = plt.subplots(figsize=(7, 4))
        names = [os.path.basename(f).replace(".family_sizes.txt", "") for f in a.family_sizes]
        scm = colour_map(names)
        for f, n in zip(a.family_sizes, names):
            d = pd.read_csv(f, sep="\t")
            d = d[d.family_size <= 30]
            ax.plot(d.family_size, d.ss_fraction, color=scm[n], label=n)
        ax.set_xlabel("Single-strand family size (reads)")
        ax.set_ylabel("Fraction of families")
        ax.set_title("UMI family-size distribution", loc="left")
        ax.legend(fontsize=8)
        save(fig, pdir, "family_sizes.png")

    if a.duplex_yield:
        fig, ax = plt.subplots(figsize=(7, 4))
        names = [os.path.basename(f).replace(".duplex_yield_metrics.txt", "") for f in a.duplex_yield]
        scm = colour_map(names)
        for f, n in zip(a.duplex_yield, names):
            d = pd.read_csv(f, sep="\t")
            ax.plot(d.read_pairs / 1e6, d.ds_duplexes / 1e3, marker="o", color=scm[n], label=n)
        ax.set_xlabel("Read pairs (millions, downsampled)")
        ax.set_ylabel("Duplex families (thousands)")
        ax.set_title("Duplex yield vs sequencing depth - where it flattens, more reads stop helping",
                     loc="left", fontsize=9)
        ax.legend(fontsize=8)
        save(fig, pdir, "duplex_yield.png")


if __name__ == "__main__":
    main()
