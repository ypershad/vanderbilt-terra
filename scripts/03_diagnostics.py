#!/usr/bin/env python3
"""
03_diagnostics.py — per-stratum QQ + Manhattan, a qc_summary.tsv, and a
meeting-ready multi-panel PDF. Reads harmonized/*.parquet (Task 1 output).

Outputs (under --fig-dir, default ../hematuria_gwas/figures):
  qq/{cohort}_{ancestry}_{phenodef}.png          300 dpi
  manhattan/{cohort}_{ancestry}_{phenodef}.png   300 dpi
  qc_summary.tsv    cohort, ancestry, pheno def, build, N pre/post QC,
                    N cases/controls, lambda_GC, lambda_1000, N loci p<5e-8
  cohort_diagnostics.pdf   one page per cohort (primary phenodef strata),
                    + two placeholder pages: "MVP" and "All of Us v9"

lambda_1000 is reported as n/a until case/control counts land in cohorts.yaml
(user decision). Everything else is complete.

Run on the VM from the repo dir, AFTER 02_harmonize.py:
  python3 scripts/03_diagnostics.py \
      --harm-dir ../hematuria_gwas/harmonized \
      --config config/cohorts.yaml \
      --data-dir ../hematuria_gwas/biovu_results \
      --fig-dir ../hematuria_gwas/figures \
      --workers 12
"""
import argparse
import os
import sys
import traceback
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")                       # headless VM
import matplotlib.pyplot as plt             # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
import polars as pl                         # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_sumstats as L                    # noqa: E402
import lib_gwasplot as G                    # noqa: E402

COHORT_ORDER = ["biovu_agd250k", "mgi", "ukb", "finngen", "taiwan", "koges"]
PLACEHOLDERS = ["MVP", "All of Us v9"]


def _sid(job):
    return f"{job['cohort']}_{job['ancestry']}_{job['phenodef']}"


def diagnose_one(task):
    """Worker: render QQ + Manhattan PNGs for one stratum, return a summary row."""
    job, fig_dir, n_raw = task["job"], task["fig_dir"], task["n_raw"]
    sid = _sid(job)
    try:
        df = pl.read_parquet(os.path.join(task["harm_dir"], f"{sid}.parquet"))
        chrom = df["CHR"].to_numpy()
        pos = df["POS"].to_numpy()
        lp = df["LOG10P"].to_numpy()
        rsid = df["RSID"].to_numpy()
        beta = df["BETA"].to_numpy()
        se = df["SE"].to_numpy()

        lam = G.lambda_gc(beta, se)
        lam1k = G.lambda_1000(lam, job.get("n_cases"), job.get("n_controls"))
        leads = G.clump_lead_loci(chrom, pos, lp, rsid)

        # QQ
        figq, axq = plt.subplots(figsize=(4.2, 4.2))
        G.render_qq(axq, lp, lam, lam1k, title=sid.replace("_", " "))
        figq.tight_layout()
        qq_path = os.path.join(fig_dir, "qq", f"{sid}.png")
        figq.savefig(qq_path, dpi=300)
        plt.close(figq)

        # Manhattan
        figm, axm = plt.subplots(figsize=(10, 3.6))
        G.render_manhattan(axm, chrom, pos, lp, rsid, leads, title=sid.replace("_", " "))
        figm.tight_layout()
        mh_path = os.path.join(fig_dir, "manhattan", f"{sid}.png")
        figm.savefig(mh_path, dpi=300)
        plt.close(figm)

        return dict(ok=True, sid=sid, row=dict(
            cohort=job["cohort"], ancestry=job["ancestry"], phenodef=job["phenodef"],
            build=job.get("build"), n_pre_qc=n_raw, n_post_qc=df.height,
            n_cases=job.get("n_cases") if job.get("n_cases") else "NA",
            n_controls=job.get("n_controls") if job.get("n_controls") else "NA",
            lambda_gc=round(lam, 4), lambda_1000=("NA" if lam1k != lam1k else round(lam1k, 4)),
            n_loci_gwsig=len(leads),
            rounded=job.get("rounded_precision", False),
        ), qq=qq_path, mh=mh_path)
    except Exception:
        return dict(ok=False, sid=sid, error=traceback.format_exc())


def build_pdf(pdf_path, jobs_by_cohort, primary_pheno, png_index):
    """One page per cohort (primary-phenodef strata) + placeholder pages."""
    with PdfPages(pdf_path) as pdf:
        for cohort in COHORT_ORDER:
            jobs = jobs_by_cohort.get(cohort, [])
            # primary strata: phenodef == primary, else all (external cohorts differ)
            prim = [j for j in jobs if j["phenodef"] == primary_pheno] or jobs
            prim = [j for j in prim if _sid(j) in png_index]
            if not prim:
                continue
            nrows = len(prim)
            fig = plt.figure(figsize=(11, 3.4 * nrows + 0.4))
            fig.suptitle(cohort, fontsize=13, y=0.995)
            for i, j in enumerate(prim):
                sid = _sid(j)
                axq = fig.add_subplot(nrows, 2, 2 * i + 1)
                axm = fig.add_subplot(nrows, 2, 2 * i + 2)
                axq.imshow(plt.imread(png_index[sid]["qq"])); axq.axis("off")
                axm.imshow(plt.imread(png_index[sid]["mh"])); axm.axis("off")
                axq.set_title(sid.replace("_", " "), fontsize=8)
            fig.tight_layout(rect=[0, 0, 1, 0.98])
            pdf.savefig(fig)
            plt.close(fig)

        # placeholder pages — user pastes MVP / AoU-v9 figures manually
        for name in PLACEHOLDERS:
            fig = plt.figure(figsize=(11, 8.5))
            fig.text(0.5, 0.5, name, ha="center", va="center", fontsize=40,
                     color="#BBBBBB")
            fig.text(0.5, 0.42, "(placeholder — paste figures here)", ha="center",
                     va="center", fontsize=12, color="#BBBBBB")
            pdf.savefig(fig)
            plt.close(fig)
    print(f"Wrote {pdf_path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--harm-dir", default="../hematuria_gwas/harmonized")
    ap.add_argument("--config", default="config/cohorts.yaml")
    ap.add_argument("--data-dir", default="../hematuria_gwas/biovu_results")
    ap.add_argument("--fig-dir", default="../hematuria_gwas/figures")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    for sub in ("qq", "manhattan"):
        os.makedirs(os.path.join(args.fig_dir, sub), exist_ok=True)

    jobs = L.resolve_configs(args.config, args.data_dir)

    # pre-QC N per stratum from harmonize_summary.tsv (fall back to NA)
    n_raw_map = {}
    hs = os.path.join(args.harm_dir, "harmonize_summary.tsv")
    if os.path.exists(hs):
        for r in pl.read_csv(hs, separator="\t").to_dicts():
            n_raw_map[f"{r['cohort']}_{r['ancestry']}_{r['phenodef']}"] = r["n_raw"]

    # only diagnose strata whose parquet exists
    tasks = []
    for j in jobs:
        p = os.path.join(args.harm_dir, f"{_sid(j)}.parquet")
        if os.path.exists(p):
            tasks.append(dict(job=j, fig_dir=args.fig_dir, harm_dir=args.harm_dir,
                              n_raw=n_raw_map.get(_sid(j), "NA")))
        else:
            print(f"skip (no parquet): {_sid(j)}")

    print(f"Diagnosing {len(tasks)} strata with {args.workers} workers ...")
    with Pool(processes=args.workers, maxtasksperchild=1) as pool:
        results = pool.map(diagnose_one, tasks)

    rows, png_index, failed = [], {}, []
    for r in results:
        if not r["ok"]:
            failed.append(r["sid"])
            print(f"\n!! FAILED {r['sid']}:\n{r['error']}")
            continue
        rows.append(r["row"])
        png_index[r["sid"]] = dict(qq=r["qq"], mh=r["mh"])
        rr = r["row"]
        print(f"  {r['sid']:<28} lambda_GC={rr['lambda_gc']:<7} "
              f"loci(5e-8)={rr['n_loci_gwsig']:<4} post_qc={rr['n_post_qc']:,}")

    if rows:
        pl.DataFrame(rows).write_csv(os.path.join(args.fig_dir, "qc_summary.tsv"),
                                     separator="\t")
        print(f"\nWrote {args.fig_dir}/qc_summary.tsv")

    # read primary phenodef from config for PDF grouping
    import yaml
    with open(args.config) as fh:
        primary_pheno = yaml.safe_load(fh)["meta_groups"].get("primary_phenodef", "593")
    jobs_by_cohort = {}
    for j in jobs:
        jobs_by_cohort.setdefault(j["cohort"], []).append(j)
    build_pdf(os.path.join(args.fig_dir, "cohort_diagnostics.pdf"),
              jobs_by_cohort, primary_pheno, png_index)

    print(f"\nDone. {len(rows)}/{len(tasks)} strata. Failed: {failed or 'none'}")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
