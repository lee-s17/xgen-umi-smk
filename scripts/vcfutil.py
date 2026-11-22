"""Small shared helpers: VCF reading, pileup AD parsing, exact tail probabilities."""
import gzip
import math


def vopen(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def info_get(info, tag):
    for kv in info.split(";"):
        if kv.startswith(tag + "="):
            return kv.split("=", 1)[1].split(",")[0]
    return None


def read_bed_set(path):
    """BED -> {chrom: sorted [(start0, end)]} for simple membership tests (gz or plain)."""
    import bisect  # noqa: F401
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    out = {}
    with (gzip.open(path, "rt") if gz else open(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c, s, e = line.split("\t")[:3]
            out.setdefault(c, []).append((int(s), int(e)))
    return {c: sorted(v) for c, v in out.items()}


def in_bed_set(bed, chrom, pos1):
    import bisect
    iv = bed.get(chrom)
    if not iv:
        return False
    i = bisect.bisect_right(iv, (pos1 - 1, float("inf"))) - 1
    return i >= 0 and iv[i][0] <= pos1 - 1 < iv[i][1]


def read_excluded(path):
    """Keys from a truth_in_normals TSV (chrom, pos, ref, alt, ...); empty if no path."""
    out = set()
    if path:
        with open(path) as fh:
            next(fh, None)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                out.add((f[0], int(f[1]), f[2], f[3]))
    return out


def read_truth(path, af_tag="", default_vaf=None):
    """{(chrom,pos,ref,alt): truth VAF or None}"""
    out = {}
    with vopen(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            alt = f[4]
            if alt in (".", "*") or alt.startswith("<"):
                continue
            af = None
            if af_tag:
                v = info_get(f[7], af_tag)
                af = float(v) if v not in (None, ".") else None
            if af is None:
                af = default_vaf
            out[(f[0], int(f[1]), f[3], alt)] = af
    return out


def read_pileup_ad(path):
    """Normalised bcftools mpileup (-a AD,DP; after norm -m-any) ->
    {(chrom,pos,ref,alt): [(alt_count, depth), ...] per sample in column order}"""
    out = {}
    with vopen(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            alt = f[4]
            if alt in (".", "*") or alt.startswith("<"):
                continue
            fmt = f[8].split(":")
            iad, idp = fmt.index("AD"), fmt.index("DP")
            per = []
            for col in f[9:]:
                v = col.split(":")
                ad = v[iad].split(",")
                a = int(ad[1]) if len(ad) > 1 and ad[1] != "." else 0
                d = int(v[idp]) if v[idp] != "." else 0
                per.append((a, d))
            out[(f[0], int(f[1]), f[3], alt)] = per
    return out


def binom_sf(k, n, p):
    """P(X >= k) for X ~ Binomial(n, p), via log-space summation."""
    if k <= 0:
        return 1.0
    if k > n or p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    lp, lq = math.log(p), math.log1p(-p)
    tot = 0.0
    for i in range(k, n + 1):
        tot += math.exp(math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) + i * lp + (n - i) * lq)
    return min(1.0, tot)


def poisson_sf(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    if lam <= 0:
        return 0.0
    cdf = sum(math.exp(-lam + i * math.log(lam) - math.lgamma(i + 1)) for i in range(k))
    return max(0.0, 1.0 - cdf)
