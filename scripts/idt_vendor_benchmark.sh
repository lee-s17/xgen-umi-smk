#!/usr/bin/env bash
# Lift over (if needed), normalise and benchmark the SEQC2 IDT site VCFs
# against the same truth set and panel as the pipeline.
#
# usage: bash scripts/idt_vendor_benchmark.sh [vcf_dir] [out_dir]
#   defaults: resources/idt2  results/idt
# env overrides: REF, PANEL_BED, NEG_BED, CHAIN, TRUTH, EXCL
set -euo pipefail

IN=${1:-resources/idt2}
OUT=${2:-results/idt}
cfg() { python -c "import yaml,sys;c=yaml.safe_load(open('config/config.yaml'));print(eval(sys.argv[1]))" "$1"; }

REF=${REF:-$(cfg 'c["ref"]')}
PANEL_BED=${PANEL_BED:-results/ref/panel.bed}
NEG_BED=${NEG_BED:-$(cfg 'c.get("negative_bed","")')}
TRUTH=${TRUTH:-results/truth/seqc2_pool.norm.vcf.gz}
EXCL=${EXCL:-results/truth/seqc2_pool.in_normals.tsv}
CHAIN=${CHAIN:-resources/liftover/hg19ToHg38.over.chain.gz}
BINS=$(cfg '",".join(map(str,c["benchmark"]["vaf_bins"]))')
HG38_CHR1=248956422

mkdir -p "$OUT/tmp" "$OUT/bench" "$(dirname "$CHAIN")"
for f in "$REF" "$PANEL_BED" "$TRUTH"; do
  [[ -e $f ]] || { echo "missing: $f" >&2; exit 1; }
done
NEG_ARG=""; [[ -n $NEG_BED && -e $NEG_BED ]] && NEG_ARG="--negative-bed $NEG_BED"
EXCL_ARG=""; [[ -e $EXCL ]] && EXCL_ARG="--exclude $EXCL"

# expected VAF scale from the sample name (fraction of Sample A)
scale_for() {
  case $1 in
    SampleA*)          echo 1 ;;
    SampleB*)          echo 0 ;;     # negative control: FPs only
    SampleC*)          echo 0.2 ;;
    SampleD*)          echo 0.2 ;;   # 1:4
    SampleE*)          echo 0.04 ;;  # 1:24 (Ef, EfIS and Ep alike)
    *)                 echo "" ;;
  esac
}

# chr1 -> chr1 rename map for VCFs using 1..22,X,Y,MT
MAP=$OUT/tmp/num2chr.txt
{ for c in $(seq 1 22) X Y; do printf '%s\tchr%s\n' "$c" "$c"; done; printf 'MT\tchrM\n'; } > "$MAP"

LOG=$OUT/liftover_counts.tsv
printf 'sample\tbuild\tn_in\tn_lifted\tn_rejected\tn_pass_norm\n' > "$LOG"

shopt -s nullglob
for vcf in "$IN"/*.vcf.gz "$IN"/*.vcf; do
  s=$(basename "$vcf"); s=${s%.gz}; s=${s%.vcf}
  scale=$(scale_for "$s")
  [[ -n $scale ]] || { echo "skip $s: unknown sample type" >&2; continue; }
  echo "== $s (vaf_scale $scale)"

  work=$OUT/tmp/$s.in.vcf.gz
  bcftools view "$vcf" -Oz -o "$work"
  # add chr prefix if needed
  if ! bcftools view -h "$work" | grep -q '^##contig=<ID=chr'; then
    if bcftools view -H "$work" | head -1 | cut -f1 | grep -qv '^chr'; then
      bcftools annotate --rename-chrs "$MAP" "$work" -Oz -o "$work.tmp" && mv "$work.tmp" "$work"
    fi
  fi
  n_in=$(bcftools view -H "$work" | wc -l)

  len=$(bcftools view -h "$work" | grep -m1 -oP '##contig=<ID=chr1,.*?length=\K[0-9]+' || true)
  if [[ $len == "$HG38_CHR1" ]]; then
    build=hg38; lifted=$work; n_rej=0
  else
    build=${len:+hg19}; build=${build:-unknown_assumed_hg19}
    [[ -e $CHAIN ]] || wget -q -O "$CHAIN" https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
    lifted=$OUT/tmp/$s.hg38.vcf.gz
    gatk --java-options "-Xmx8g" LiftoverVcf -I "$work" -O "$lifted" -C "$CHAIN" -R "$REF" \
      --REJECT "$OUT/tmp/$s.rejected.vcf.gz" --RECOVER_SWAPPED_REF_ALT true \
      --WARN_ON_MISSING_CONTIG true --VERBOSITY WARNING 2> "$OUT/tmp/$s.liftover.log"
    n_rej=$(bcftools view -H "$OUT/tmp/$s.rejected.vcf.gz" | wc -l)
  fi
  n_lift=$(bcftools view -H "$lifted" | wc -l)

  final=$OUT/$s.idt.hg38.vcf.gz
  bcftools reheader --fai "$REF.fai" "$lifted" \
    | bcftools sort -Ou - \
    | bcftools norm -a -m-any -f "$REF" -c w -Ou - \
    | bcftools view -f PASS,. -Oz -o "$final"
  bcftools index -tf "$final"
  n_pass=$(bcftools view -H "$final" | wc -l)
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$s" "$build" "$n_in" "$n_lift" "$n_rej" "$n_pass" >> "$LOG"

  python scripts/benchmark.py --calls "$final" --truth "$TRUTH" --bed "$PANEL_BED" $NEG_ARG \
    --vaf-scale "$scale" --truth-af-tag VAF --default-truth-vaf 0 --bins "$BINS" \
    --sample "$s" --level idt_vendor $EXCL_ARG --out "$OUT/bench/$s.tsv"
done

# one long table: all samples, all metrics
first=1
for t in "$OUT"/bench/*.tsv; do
  if (( first )); then cat "$t"; first=0; else tail -n +2 "$t"; fi
done > "$OUT/idt_vendor_benchmark.tsv"
echo
column -t "$LOG"
echo "benchmark table: $OUT/idt_vendor_benchmark.tsv"
