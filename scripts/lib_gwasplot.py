"""
lib_gwasplot.py — GWAS diagnostic-plot helpers (QQ, Manhattan) + genetics utils.

Shared by 03_diagnostics.py (per-stratum) and the Task-3 post-meta plots. Static
matplotlib (Agg) at 300 dpi for slide/screen-share. Genre conventions of a genetics
audience govern here; the CVD-safe accent (significant hits vs chromosome bands) was
validated with the dataviz palette checker (ΔE 20.7 deutan, well above the ≥8 target).

Design points:
  - Works off the harmonized LOG10P column (underflow-safe) — no inf at GWAS sig.
  - lambda_GC from Wald chisq=(BETA/SE)^2 (uniform across all cohorts; robust at the
    median, so p-underflow is irrelevant).
  - Null thinning: keep every point with -log10p above a threshold, subsample the
    bulk. Files stay small without distorting the tail.
"""
import numpy as np

# --- palette (validated: accent vs bands CVD ΔE 20.7; bands = recessive parity) ---
BAND_DARK = "#2F4B82"
BAND_LIGHT = "#6E93C6"
ACCENT = "#E1622F"       # genome-wide significant hits
GW_LINE = "#C1272D"      # p = 5e-8
SUG_LINE = "#8A8A8A"     # p = 1e-5
INK = "#222222"

GW_LOG10P = -np.log10(5e-8)     # 7.30103
SUG_LOG10P = -np.log10(1e-5)    # 5.0
QCHISQ_MED_1DF = 0.4549364      # qchisq(0.5, df=1)

CHR_ORDER = [str(i) for i in range(1, 23)] + ["X"]


# --------------------------------------------------------------- genetics utils
def lambda_gc(beta, se):
    """Genomic inflation from Wald chisq. NaN-safe; returns float."""
    beta = np.asarray(beta, float)
    se = np.asarray(se, float)
    ok = np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if ok.sum() == 0:
        return float("nan")
    chisq = (beta[ok] / se[ok]) ** 2
    return float(np.median(chisq) / QCHISQ_MED_1DF)


def lambda_1000(lam, n_cases, n_controls):
    """Case/control-size-corrected inflation, scaled to 1000/1000.

    lambda_1000 = 1 + (lambda-1) * (1/Ncase + 1/Nctrl) / (1/1000 + 1/1000).
    Returns NaN if counts are missing (deferred until counts are supplied).
    """
    if not n_cases or not n_controls or n_cases <= 0 or n_controls <= 0:
        return float("nan")
    return 1.0 + (lam - 1.0) * (1.0 / n_cases + 1.0 / n_controls) / (2.0 / 1000.0)


def clump_lead_loci(chrom, pos, log10p, rsid=None, window=1_000_000,
                    p_thresh=5e-8):
    """Greedy 1 Mb clumping on p. Returns list of lead dicts, most-sig first.

    Only variants with p < p_thresh are eligible leads. For each, remove all
    variants within +/- window on the same chromosome, repeat.
    """
    thr = -np.log10(p_thresh)
    chrom = np.asarray(chrom)
    pos = np.asarray(pos, float)
    lp = np.asarray(log10p, float)
    rsid = np.asarray(rsid) if rsid is not None else np.array([None] * len(pos))

    sig = np.where(np.isfinite(lp) & (lp > thr))[0]
    if sig.size == 0:
        return []
    order = sig[np.argsort(-lp[sig])]        # descending significance
    taken = np.zeros(len(pos), bool)
    leads = []
    for i in order:
        if taken[i]:
            continue
        leads.append(dict(chr=str(chrom[i]), pos=int(pos[i]),
                          log10p=float(lp[i]),
                          rsid=(None if rsid[i] in (None, "", "nan") else str(rsid[i]))))
        same = (chrom == chrom[i]) & (np.abs(pos - pos[i]) <= window)
        taken |= same
    return leads


def _chrom_layout(chrom, pos, gap_frac=0.02):
    """Cumulative x-position per chromosome + tick centers. Returns (x, ticks, labels)."""
    chrom = np.asarray(chrom).astype(str)
    pos = np.asarray(pos, float)
    present = [c for c in CHR_ORDER if (chrom == c).any()]
    x = np.zeros(len(pos))
    offset = 0.0
    ticks, labels = [], []
    span_total = max(1.0, sum(pos[chrom == c].max() for c in present))
    gap = span_total * gap_frac / max(1, len(present))
    for c in present:
        m = chrom == c
        cmax = pos[m].max()
        x[m] = offset + pos[m]
        ticks.append(offset + cmax / 2.0)
        labels.append(c)
        offset += cmax + gap
    return x, ticks, labels, present


def _thin(log10p, keep_above=4.0, n_bulk=200_000, seed_stride=None):
    """Index mask: keep all points above `keep_above`, subsample the rest.

    Deterministic stride subsample (no RNG — reproducible for methods).
    """
    lp = np.asarray(log10p, float)
    keep = np.where(np.isfinite(lp) & (lp > keep_above))[0]
    bulk = np.where(np.isfinite(lp) & (lp <= keep_above))[0]
    if bulk.size > n_bulk:
        stride = seed_stride or max(1, bulk.size // n_bulk)
        bulk = bulk[::stride]
    return np.sort(np.concatenate([keep, bulk]))


# --------------------------------------------------------------------- renders
def render_qq(ax, log10p, lam_gc, lam_1000=float("nan"), title=""):
    """QQ of observed vs expected -log10 p, with lambda printed. Null thinned."""
    lp = np.asarray(log10p, float)
    lp = lp[np.isfinite(lp)]
    obs = np.sort(lp)[::-1]                       # descending: obs[0] = smallest p
    M = obs.size
    exp = -np.log10((np.arange(1, M + 1) - 0.5) / M)

    sel = _thin(obs, keep_above=4.0, n_bulk=15_000)   # thin the dense null
    lim = float(max(exp.max(), obs.max())) * 1.05

    ax.plot([0, lim], [0, lim], color=SUG_LINE, lw=1, zorder=1)
    ax.scatter(exp[sel], obs[sel], s=6, color=BAND_DARK, edgecolors="none", zorder=2)
    ax.axhline(GW_LOG10P, color=GW_LINE, lw=0.8, ls="--", zorder=1)

    lam_txt = f"$\\lambda_{{GC}}$ = {lam_gc:.3f}"
    if np.isfinite(lam_1000):
        lam_txt += f"\n$\\lambda_{{1000}}$ = {lam_1000:.3f}"
    else:
        lam_txt += "\n$\\lambda_{1000}$ = n/a"
    ax.text(0.05, 0.95, lam_txt, transform=ax.transAxes, va="top", ha="left",
            fontsize=9, color=INK,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=SUG_LINE, lw=0.6))
    ax.set_xlabel("Expected $-\\log_{10}(p)$", fontsize=9)
    ax.set_ylabel("Observed $-\\log_{10}(p)$", fontsize=9)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    ax.tick_params(labelsize=8)


def render_manhattan(ax, chrom, pos, log10p, rsid=None, leads=None,
                     title="", max_labels=20):
    """Manhattan with 5e-8 / 1e-5 lines and labeled lead SNPs. Null thinned."""
    chrom = np.asarray(chrom).astype(str)
    pos = np.asarray(pos, float)
    lp = np.asarray(log10p, float)

    sel = _thin(lp, keep_above=3.0, n_bulk=300_000)
    x_all, ticks, labels, present = _chrom_layout(chrom, pos)

    # plot thinned, colored by chromosome parity (recessive bands)
    parity = {c: i % 2 for i, c in enumerate(present)}
    cs = chrom[sel]
    colors = np.where(np.vectorize(parity.get)(cs) == 0, BAND_DARK, BAND_LIGHT)
    ax.scatter(x_all[sel], lp[sel], s=4, c=colors, edgecolors="none", zorder=2)

    # significant hits in accent on top
    sig = sel[lp[sel] > GW_LOG10P]
    if sig.size:
        ax.scatter(x_all[sig], lp[sig], s=9, color=ACCENT, edgecolors="none", zorder=3)

    ax.axhline(GW_LOG10P, color=GW_LINE, lw=0.8, ls="--", zorder=1)
    ax.axhline(SUG_LOG10P, color=SUG_LINE, lw=0.8, ls=":", zorder=1)

    # label lead SNPs (already clumped); map lead back to its cumulative x
    if leads:
        # build a chrom->offset lookup (same layout as _chrom_layout) for lead x
        off = {}
        acc = 0.0
        span_total = max(1.0, sum(pos[chrom == c].max() for c in present))
        gap = span_total * 0.02 / max(1, len(present))
        for c in present:
            off[c] = acc
            acc += pos[chrom == c].max() + gap
        for ld in leads[:max_labels]:
            if ld["chr"] not in off:
                continue
            lx = off[ld["chr"]] + ld["pos"]
            lab = ld["rsid"] or f"{ld['chr']}:{ld['pos']}"
            ax.annotate(lab, (lx, ld["log10p"]), fontsize=6.5, color=INK,
                        xytext=(0, 4), textcoords="offset points",
                        ha="center", rotation=45)

    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_xlim(x_all.min(), x_all.max())
    ax.set_ylim(0, max(GW_LOG10P + 1, float(np.nanmax(lp)) * 1.05))
    ax.set_xlabel("Chromosome", fontsize=9)
    ax.set_ylabel("$-\\log_{10}(p)$", fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, color=INK)
    ax.tick_params(axis="y", labelsize=8)
