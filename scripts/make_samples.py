#!/usr/bin/env python3
"""Build config/samples.tsv from the SEQC2 IDT experiment list (+ optional ENA run table).

Library names look like SampleEfIS_IDT2_ST04_25ng_LIB1:
  B/D/E = sample, f = fragmented, IS = Accugenomics spike-in, p = synthetic plasma (Ep)
  ST04-06 = site, 10/25/50ng = input, LIB1-4 = technical replicate

Get the run table on a machine with internet:
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=SRP296025&result=read_run&fields=run_accession,experiment_alias,read_count,base_count&format=tsv&limit=0" > resources/metadata/SRP296025_runs.tsv

Examples:
  python scripts/make_samples.py --runs resources/metadata/SRP296025_runs.tsv --fill config/samples.tsv  # fill SRRs only
  python scripts/make_samples.py --runs resources/metadata/SRP296025_runs.tsv --reps 2 --titration  # default subset (8 libs)
  python scripts/make_samples.py --runs resources/metadata/SRP296025_runs.tsv                      # core: ST04, 25 ng, 4 reps
  python scripts/make_samples.py --runs ... --site ST05 --titration --plasma                     # core + extras
  python scripts/make_samples.py --runs ... --site all                                            # all 3 sites
"""
import argparse
import re
import sys

import pandas as pd

PAT = re.compile(r"^Sample(?P<base>[BDE])(?P<kind>fIS|p)_IDT2_(?P<site>ST\d+)_(?P<ng>\d+)ng_(?P<lib>LIB\d+)$")
# fraction of Sample A in the mix (Gong et al. 2022: A:B = 1:4 -> D, 1:24 -> E)
VAF_SCALE = {"D": "0.2", "E": "0.04", "B": ""}


def parse(alias):
    m = PAT.match(alias)
    if not m:
        sys.exit(f"cannot parse library name: {alias}")
    d = m.groupdict()
    d["ng"] = int(d["ng"])
    d["group"] = "Ep" if d["kind"] == "p" else f"{d['base']}fIS"
    return d


ap = argparse.ArgumentParser()
ap.add_argument("--experiments", default="resources/metadata/SRP296025_IDT_experiments.tsv")
ap.add_argument("--runs", help="ENA filereport TSV (run_accession, experiment_alias)")
ap.add_argument("--site", default="ST04", help="ST04, ST05, ST06 or all")
ap.add_argument("--titration", action="store_true", help="add EfIS 10 and 50 ng")
ap.add_argument("--plasma", action="store_true", help="add Ep 25 ng")
ap.add_argument("--reps", type=int, default=4, help="replicates per 25 ng group (LIB1..N)")
ap.add_argument("--titration-reps", type=int, default=1, help="replicates per 10/50 ng input")
ap.add_argument("--truth", default="seqc2_pool")
ap.add_argument("-o", "--out", default="config/samples.tsv")
ap.add_argument("--fill", help="only fill the srr column of this existing samples.tsv (keeps comments and rows)")
a = ap.parse_args()

if a.fill:
    if not a.runs:
        sys.exit("--fill needs --runs")
    runs = pd.read_csv(a.runs, sep="\t", dtype=str)
    srr = runs.sort_values("run_accession").groupby("experiment_alias").run_accession.agg(";".join).to_dict()
    lines, n, miss = [], 0, []
    for line in open(a.fill):
        f = line.rstrip("\n").split("\t")
        if not line.startswith("#") and f[0] != "sample" and len(f) > 1:
            if f[0] in srr:
                f[1] = srr[f[0]]
                n += 1
            else:
                miss.append(f[0])
            line = "\t".join(f) + "\n"
        lines.append(line)
    open(a.fill, "w").writelines(lines)
    print(f"filled {n} rows in {a.fill}; missing: {miss or 'none'}", file=sys.stderr)
    sys.exit(0)

exp = pd.read_csv(a.experiments, sep="\t", dtype=str).dropna(subset=["experiment_alias"])
exp = exp[exp.experiment_alias.str.contains("IDT2")]
meta = pd.DataFrame([parse(x) for x in exp.experiment_alias], index=exp.experiment_alias)

keep = (meta.ng == 25) & meta.group.isin(["BfIS", "DfIS", "EfIS"])
if a.titration:
    keep |= (meta.group == "EfIS") & meta.ng.isin([10, 50])
if a.plasma:
    keep |= (meta.group == "Ep") & (meta.ng == 25)
if a.site != "all":
    keep &= meta.site == a.site
libn = meta.lib.str[3:].astype(int)
keep &= ((meta.ng == 25) & (libn <= a.reps)) | ((meta.ng != 25) & (libn <= a.titration_reps))
sel = meta[keep].copy()
if sel.empty:
    sys.exit("no libraries selected; check --site")

# one experiment can have several runs: join with ';' (fetch_sra concatenates them)
srr = {}
if a.runs:
    runs = pd.read_csv(a.runs, sep="\t", dtype=str)
    runs = runs[runs.experiment_alias.isin(sel.index)]
    srr = runs.sort_values("run_accession").groupby("experiment_alias").run_accession.agg(";".join).to_dict()
    missing = [x for x in sel.index if x not in srr]
    if missing:
        print(f"WARNING: no run found for {len(missing)} libraries, e.g. {missing[:3]}", file=sys.stderr)

order = {"DfIS": 0, "EfIS": 1, "Ep": 2, "BfIS": 3}
sel["o"] = sel.group.map(order)
sel = sel.sort_values(["site", "o", "ng", "lib"])

rows = []
for alias, r in sel.iterrows():
    normal = r.base == "B"
    rows.append({
        "sample": alias,
        "srr": srr.get(alias, "SRR_FILL_ME"),
        "fq1": "", "fq2": "",
        "truth": "" if normal else a.truth,
        "vaf_scale": VAF_SCALE[r.base],
        "role": "normal" if normal else "test",
    })

with open(a.out, "w") as fh:
    fh.write("# generated by scripts/make_samples.py from the SRP296025 IDT experiment list\n")
    fh.write("# D = Lbx-high (A:B 1:4, vaf_scale 0.2); E/Ep = Lbx-low (1:24, 0.04); B = background -> normals\n")
    pd.DataFrame(rows).to_csv(fh, sep="\t", index=False)

print(sel.groupby(["site", "group", "ng"]).size().to_string(), file=sys.stderr)
print(f"wrote {len(rows)} libraries to {a.out}", file=sys.stderr)
