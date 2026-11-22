#!/usr/bin/env bash
# Build a BED of spike-in / Sample-B-only sites to mask from benchmarking.
# Any call at a known-negative panel position with VAF >= MIN_AF in a Sample B
# library (no Sample A material) is either an Accugenomics spike-in or a
# truth-set gap, so it is masked for every pipeline alike.
#
# usage: bash scripts/make_spikein_mask.sh out.bed vcf1 [vcf2 ...]
# env:   PANEL_BED NEG_BED MIN_AF (default 0.05) PAD (default 1)
set -euo pipefail
OUT=$1; shift
PANEL_BED=${PANEL_BED:-results/ref/panel.bed}
NEG_BED=${NEG_BED:-resources/truth/KnownNegatives_hg38.bed.gz}
MIN_AF=${MIN_AF:-0.05}
PAD=${PAD:-1}

tmp=$(mktemp)
for v in "$@"; do
  # VAF from FORMAT/AF if declared, else INFO/AF
  if bcftools view -h "$v" | grep -q '^##FORMAT=<ID=AF,'; then fmt='[%AF]'
  elif bcftools view -h "$v" | grep -q '^##INFO=<ID=AF,'; then fmt='%INFO/AF'
  else echo "skip $v: no AF field" >&2; continue; fi
  # all records, any FILTER: our normals keep spike-ins as FILTER=PoN
  bcftools view -T "$PANEL_BED" "$v" -Ou \
    | bcftools view -T "$NEG_BED" -Ou \
    | bcftools query -f "%CHROM\t%POS\t$fmt\n" \
    | awk -v m="$MIN_AF" -v p="$PAD" 'BEGIN{OFS="\t"}
        { split($3,a,","); if (a[1]+0 >= m) print $1,$2-1-p,$2+p }' >> "$tmp"
  echo "$v: $(wc -l < "$tmp") sites so far" >&2
done
# benchmark.py merges overlapping intervals itself, so no bedtools needed
sort -k1,1 -k2,2n -u "$tmp" | awk 'BEGIN{OFS="\t"} {if ($2<0) $2=0; print}' > "$OUT"
rm -f "$tmp"
echo "mask: $(wc -l < "$OUT") intervals -> $OUT" >&2
