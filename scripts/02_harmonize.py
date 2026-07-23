#!/usr/bin/env python3
"""
02_harmonize.py — harmonize all 13 hematuria sumstat files to a common GRCh38
schema, in parallel (one worker per file), with every filter count logged.

Pipeline per file:
  read -> standardize(column_map) -> LOG10P (underflow-safe) -> key38
       -> liftover 37->38 (only cohorts flagged by 01_build_check) -> QC
       -> write harmonized/{cohort}_{ancestry}_{phenodef}.parquet

Outputs:
  harmonized/*.parquet                 one per stratum (16 total: 8 BioVU + 5 ext)
  harmonized/harmonize_log.tsv         long-format filter counts (methods-ready)
  harmonized/harmonize_summary.tsv     one row per stratum (pre/post N, X count, ...)

Run on the VM from the repo dir, AFTER 01_build_check.py:
  python3 scripts/02_harmonize.py \
      --data-dir ../hematuria_gwas/biovu_results \
      --config config/cohorts.yaml \
      --build-check config/build_check.tsv \
      --chain ../hematuria_gwas/ref/hg19ToHg38.over.chain.gz \
      --out-dir ../hematuria_gwas/harmonized \
      --workers 12

MEMORY: MGI decompresses to ~15-20 GB; 12 workers peak ~70-75 GB (fits 104 GB with
headroom). maxtasksperchild=1 frees each worker's memory between files. If you OOM,
drop --workers to 8.
"""
import argparse
import os
import shutil
import sys
import traceback
from multiprocessing import Pool

# Keep each worker's polars modest so 12 workers don't oversubscribe 16 cores.
os.environ.setdefault("POLARS_MAX_THREADS", "2")

import polars as pl   # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lib_sumstats as L   # noqa: E402


def _stratum_id(job):
    return f"{job['cohort']}_{job['ancestry']}_{job['phenodef']}"


def harmonize_one(job):
    """Worker: harmonize a single file. Returns a result dict (picklable)."""
    sid = _stratum_id(job)
    log = []
    try:
        df = L.read_raw(job["path"])
        raw_n = df.height
        df = L.standardize(df, job, job["cohort"], job["ancestry"], job["phenodef"])
        df = L.add_log10p(df)
        df = L.add_key38(df)

        if job.get("liftover") is True:
            df = L.liftover_pos(df, job["chain_path"], log)
            df = L.add_key38(df)   # POS changed -> refresh key (alleles unchanged)

        df = L.qc(df, job["qc"], log)

        auto, xdf = L.split_autosome_x(df)
        n_x = xdf.height
        df = df.select(L.OUT_COLS)

        out_path = os.path.join(job["out_dir"], f"{sid}.parquet")
        df.write_parquet(out_path)

        summary = dict(
            cohort=job["cohort"], ancestry=job["ancestry"], phenodef=job["phenodef"],
            build_used=("GRCh38" if not job.get("liftover") else "GRCh38(lifted)"),
            n_raw=raw_n, n_post_qc=df.height, n_autosome=auto.height, n_chrX=n_x,
            liftover=bool(job.get("liftover")),
            has_eaf=job["column_map"].get("eaf") is not None,
            has_n=job["column_map"].get("n") is not None,
            rounded=job.get("format") in ("pheweb", "pheweb_minimal"),
            n_gwsig=int((df["LOG10P"] > 7.30103).sum()),   # p < 5e-8
            out=out_path,
        )
        return dict(sid=sid, ok=True, log=log, summary=summary)
    except Exception:
        return dict(sid=sid, ok=False, error=traceback.format_exc(), log=log)


def preflight():
    missing = []
    for mod in ("numpy", "scipy", "yaml", "polars"):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    if missing:
        sys.exit(f"Missing required packages: {missing}. pip install them first.")


def apply_build_check(jobs, path):
    """Override each suspect cohort's build/liftover from 01_build_check.py output."""
    if not os.path.exists(path):
        print(f"WARNING: {path} not found -> using yaml build flags as-is. "
              f"Run 01_build_check.py first for UKB/KoGES/Taiwan.")
        return jobs
    bc = {r["cohort"]: r for r in pl.read_csv(path, separator="\t").to_dicts()}
    for j in jobs:
        r = bc.get(j["cohort"])
        if r:
            j["liftover"] = str(r["liftover"]).lower() == "true"
            j["build"] = r["verdict"]
            if r["verdict"].startswith("AMBIGUOUS"):
                print(f"!! {j['cohort']}: build AMBIGUOUS ({r['verdict']}). "
                      f"Refusing to guess — resolve before harmonizing this cohort.")
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="../hematuria_gwas/biovu_results")
    ap.add_argument("--config", default="config/cohorts.yaml")
    ap.add_argument("--build-check", default="config/build_check.tsv")
    ap.add_argument("--chain", default="../hematuria_gwas/ref/hg19ToHg38.over.chain.gz")
    ap.add_argument("--out-dir", default="../hematuria_gwas/harmonized")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()

    preflight()
    os.makedirs(args.out_dir, exist_ok=True)

    # Disk-space guard (brief: warn if < ~200 GB free where parquet is written).
    free_gb = shutil.disk_usage(args.out_dir).free / 1e9
    print(f"Free disk at {args.out_dir}: {free_gb:.0f} GB")
    if free_gb < 200:
        print(f"!! WARNING: < 200 GB free ({free_gb:.0f} GB). Parquet outputs are far "
              f"smaller than the gzipped inputs, but confirm before proceeding.")

    jobs = L.resolve_configs(args.config, args.data_dir)
    jobs = apply_build_check(jobs, args.build_check)

    need_lift = [j["cohort"] for j in jobs if j.get("liftover") is True]
    if need_lift:
        if not os.path.exists(args.chain):
            sys.exit(f"Liftover needed for {sorted(set(need_lift))} but chain file "
                     f"not found: {args.chain}. Install pyliftover + download the "
                     f"UCSC hg19ToHg38.over.chain.gz.")
        try:
            import pyliftover  # noqa: F401
        except ImportError:
            sys.exit("Liftover needed but pyliftover not installed: pip install pyliftover")

    for j in jobs:
        j["out_dir"] = args.out_dir
        j["chain_path"] = args.chain

    print(f"Harmonizing {len(jobs)} strata with {args.workers} workers "
          f"(lift: {sorted(set(need_lift)) or 'none'}) ...")

    with Pool(processes=args.workers, maxtasksperchild=1) as pool:
        results = pool.map(harmonize_one, jobs)

    # ---- aggregate logs ----
    log_rows, summ_rows, failed = [], [], []
    for r in results:
        if not r["ok"]:
            failed.append(r["sid"])
            print(f"\n!! FAILED {r['sid']}:\n{r['error']}")
            continue
        for step in r["log"]:
            log_rows.append(dict(stratum=r["sid"], **step))
        summ_rows.append(r["summary"])
        s = r["summary"]
        print(f"  {r['sid']:<28} raw={s['n_raw']:>10,} -> post_qc={s['n_post_qc']:>10,}"
              f"  X={s['n_chrX']:>7,}  gwsig={s['n_gwsig']:>4}"
              f"  {'LIFTED' if s['liftover'] else ''}")

    if log_rows:
        pl.DataFrame(log_rows).write_csv(
            os.path.join(args.out_dir, "harmonize_log.tsv"), separator="\t")
    if summ_rows:
        pl.DataFrame(summ_rows).write_csv(
            os.path.join(args.out_dir, "harmonize_summary.tsv"), separator="\t")

    print(f"\nDone. {len(summ_rows)}/{len(jobs)} strata harmonized. "
          f"Failed: {failed or 'none'}")
    print(f"Logs: {args.out_dir}/harmonize_log.tsv, harmonize_summary.tsv")
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
