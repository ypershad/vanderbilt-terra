#!/usr/bin/env python3
"""
quicklook_plots.py — FAST per-cohort QQ + Manhattan straight from the RAW files,
skipping harmonization. For an eyeball while 02_harmonize.py runs in another tab.

Why this is safe to skip harmonization for (diagnostics only, NOT the meta):
  - A per-cohort QQ/Manhattan needs only CHR, POS, and a p-value in that cohort's
    OWN coordinates. Liftover, allele orientation, EAF/INFO/N, and the common
    cross-cohort schema are all irrelevant to a single cohort's diagnostic plot.
What it still does right (the one real gotcha):
  - P-values underflow to 0 at genome-wide significance (BioVU esp.). It reuses the
    underflow-safe LOG10P logic (FinnGen mlogp -> BioVU Chisq -> P -> Wald z), so the
    strongest loci render instead of getting clipped.
What it deliberately does NOT do (vs 03_diagnostics.py):
  - No QC filtering beyond a light sanity filter (valid chr, finite SE>0, finite
    LOG10P). So lambda/shape are a QUICK LOOK, not the methods-grade numbers. Titles
    say "(quick-look)" so nobody confuses these with the QC'd figures.
  - No qc_summary.tsv, no PDF. Those come from 03_diagnostics on harmonized parquet.

Reads only the columns a plot needs (projection) -> light memory, so it can run
alongside a harmonize job. Keep --workers low (2-3) while harmonize is mid-MGI.

Run on the VM from the repo dir (separate tab):
  python3 scripts/quicklook_plots.py \
      --data-dir ../hematuria_gwas/biovu_results \
      --config config/cohorts.yaml \
      --fig-dir ../hematuria_gwas/figures_quicklook \
      --workers 3
Optionally limit cohorts:  --only ukb,finngen,taiwan
"""
import argparse
import os
import sys
import traceback
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402
import polars as pl               # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_sumstats as L          # noqa: E402
import lib_gwasplot as G          # noqa: E402


def _sid(job):
    return f"{job['cohort']}_{job['ancestry']}_{job['phenodef']}"


def _needed_sources(cm, available):
    """Source columns a plot needs, intersected with what the file actually has."""
    keys = ["chr", "pos", "beta", "se", "p", "log10p", "chisq"]
    cols = [cm[k] for k in keys if cm.get(k)]
    if "rsids" in available:
        cols.append("rsids")       # nicer lead labels; else fall back to chr:pos
    # de-dup, preserve order, keep only columns present in the header
    seen, out = set(), []
    for c in cols:
        if c in available and c not in seen:
            seen.add(c)
            out.append(c)
    return out


def quicklook_one(task):
    job, fig_dir = task["job"], task["fig_dir"]
    sid = _sid(job)
    cm = job["column_map"]
    try:
        # header-only read to learn which columns exist, then projected read
        avail = pl.read_csv(job["path"], separator="\t", n_rows=0,
                            infer_schema_length=0).columns
        cols = _needed_sources(cm, avail)
        df = pl.read_csv(job["path"], separator="\t", columns=cols,
                         infer_schema_length=0)   # projection -> light memory

        def src_or_null(key):
            name = cm.get(key)
            if name and name in df.columns:
                return pl.col(name).cast(pl.Float64, strict=False)
            return pl.lit(None, pl.Float64)

        mini = df.select(
            L._clean_chr(pl.col(cm["chr"])).alias("CHR"),
            pl.col(cm["pos"]).cast(pl.Int64, strict=False).alias("POS"),
            src_or_null("beta").alias("BETA"),
            src_or_null("se").alias("SE"),
            src_or_null("p").alias("P"),
            src_or_null("log10p").alias("_LOG10P_SRC"),
            src_or_null("chisq").alias("_CHISQ"),
            (pl.col("rsids").cast(pl.Utf8) if "rsids" in df.columns
             else pl.lit(None, pl.Utf8)).alias("RSID"),
        )
        mini = L.add_log10p(mini)

        # light sanity filter only (NOT the full QC)
        before = mini.height
        mini = mini.filter(
            pl.col("CHR").is_in(list(L.VALID_CHR))
            & (pl.col("POS") > 0)
            & pl.col("LOG10P").is_finite()
            & pl.col("SE").is_finite() & (pl.col("SE") > 0)
            & pl.col("BETA").is_finite()
        )

        chrom = mini["CHR"].to_numpy()
        pos = mini["POS"].to_numpy()
        lp = mini["LOG10P"].to_numpy()
        rsid = mini["RSID"].to_numpy()
        beta = mini["BETA"].to_numpy()
        se = mini["SE"].to_numpy()

        lam = G.lambda_gc(beta, se)
        leads = G.clump_lead_loci(chrom, pos, lp, rsid)

        figq, axq = plt.subplots(figsize=(4.2, 4.2))
        G.render_qq(axq, lp, lam, title=f"{sid.replace('_', ' ')} (quick-look)")
        figq.tight_layout()
        figq.savefig(os.path.join(fig_dir, "qq", f"{sid}.png"), dpi=200)
        plt.close(figq)

        figm, axm = plt.subplots(figsize=(10, 3.6))
        G.render_manhattan(axm, chrom, pos, lp, rsid, leads,
                           title=f"{sid.replace('_', ' ')} (quick-look)")
        figm.tight_layout()
        figm.savefig(os.path.join(fig_dir, "manhattan", f"{sid}.png"), dpi=200)
        plt.close(figm)

        return dict(ok=True, sid=sid, lam=round(lam, 4), n_loci=len(leads),
                    n_plotted=mini.height, n_dropped=before - mini.height)
    except Exception:
        return dict(ok=False, sid=sid, error=traceback.format_exc())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="../hematuria_gwas/biovu_results")
    ap.add_argument("--config", default="config/cohorts.yaml")
    ap.add_argument("--fig-dir", default="../hematuria_gwas/figures_quicklook")
    ap.add_argument("--workers", type=int, default=3,
                    help="keep low (2-3) while a harmonize job is running")
    ap.add_argument("--only", default="", help="comma-list of cohorts to limit to")
    args = ap.parse_args()

    for sub in ("qq", "manhattan"):
        os.makedirs(os.path.join(args.fig_dir, sub), exist_ok=True)

    jobs = L.resolve_configs(args.config, args.data_dir)
    if args.only:
        want = {c.strip() for c in args.only.split(",")}
        jobs = [j for j in jobs if j["cohort"] in want]
    tasks = [dict(job=j, fig_dir=args.fig_dir) for j in jobs
             if os.path.exists(j["path"])]

    print(f"Quick-look: {len(tasks)} strata, {args.workers} workers "
          f"-> {args.fig_dir}")
    print("(minimal filtering — lambda/shape are a QUICK LOOK, not methods-grade)\n")

    with Pool(processes=args.workers, maxtasksperchild=1) as pool:
        results = pool.map(quicklook_one, tasks)

    failed = []
    for r in sorted(results, key=lambda x: x["sid"]):
        if not r["ok"]:
            failed.append(r["sid"])
            print(f"!! FAILED {r['sid']}:\n{r['error']}")
            continue
        print(f"  {r['sid']:<28} lambda_GC={r['lam']:<7} loci(5e-8)={r['n_loci']:<4} "
              f"plotted={r['n_plotted']:,} (dropped {r['n_dropped']:,})")

    print(f"\nDone. {len(results) - len(failed)}/{len(tasks)} plotted. "
          f"Failed: {failed or 'none'}")
    print(f"PNGs: {args.fig_dir}/qq/, {args.fig_dir}/manhattan/")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
