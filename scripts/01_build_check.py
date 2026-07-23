#!/usr/bin/env python3
"""
01_build_check.py — empirically determine the genome build of the cohorts whose
build is uncertain (UKB, KoGES, Taiwan), instead of trusting filenames.

Method: FinnGen R7 is confirmed GRCh38 and carries rsids. For each suspect cohort
we intersect on rsID and compare POS:
    - most positions EQUAL  -> cohort is GRCh38 (no liftover)
    - most positions DIFFER -> cohort is GRCh37 (liftover needed)
We report the match fraction and a verdict per cohort, and write
config/build_check.tsv which 02_harmonize.py reads to decide liftover.

Build is a property of coordinates, not ancestry, so FinnGen (EUR) is a valid
reference for EAS cohorts. Only variants with a non-empty rsID are used.

Run on the VM from the repo dir:
    python3 scripts/01_build_check.py \
        --data-dir ../hematuria_gwas/biovu_results \
        --config config/cohorts.yaml \
        --out config/build_check.tsv
"""
import argparse
import os
import sys

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_sumstats import read_raw, resolve_configs   # noqa: E402

SUSPECTS = ["ukb", "koges", "taiwan"]
REFERENCE = "finngen"
N_SAMPLE = 200_000            # rsIDs to compare (plenty for a clear verdict)
MATCH_38_THRESHOLD = 0.90     # >= -> GRCh38 ; <= 0.10 -> GRCh37 ; between -> ambiguous


def _rsid_pos(job, limit=None):
    """Return polars df [RSID, POS] for variants with a real rsID (rs...)."""
    cm = job["column_map"]
    df = read_raw(job["path"])
    rsid_col = job.get("rsid_col") or ("rsids" if "rsids" in df.columns else None)
    if rsid_col is None:
        return pl.DataFrame({"RSID": [], "POS": []},
                            schema={"RSID": pl.Utf8, "POS": pl.Int64})
    out = (
        df.select(
            pl.col(rsid_col).cast(pl.Utf8).alias("RSID"),
            pl.col(cm["pos"]).cast(pl.Int64, strict=False).alias("POS"),
        )
        .filter(pl.col("RSID").str.starts_with("rs") & pl.col("POS").is_not_null())
        # PheWeb sometimes packs multiple rsids as 'rs1,rs2' -> take the first
        .with_columns(pl.col("RSID").str.split(",").list.first())
        .unique(subset="RSID")
    )
    if limit:
        out = out.head(limit)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="../hematuria_gwas/biovu_results")
    ap.add_argument("--config", default="config/cohorts.yaml")
    ap.add_argument("--out", default="config/build_check.tsv")
    args = ap.parse_args()

    jobs = {j["cohort"]: j for j in resolve_configs(args.config, args.data_dir)}

    print(f"Loading reference {REFERENCE} (GRCh38) rsID->POS ...")
    ref = _rsid_pos(jobs[REFERENCE]).rename({"POS": "POS_REF"})
    print(f"  reference rsIDs: {ref.height:,}")

    rows = []
    for name in SUSPECTS:
        print(f"\nChecking {name} ...")
        sus = _rsid_pos(jobs[name], limit=N_SAMPLE).rename({"POS": "POS_SUS"})
        j = sus.join(ref, on="RSID", how="inner")
        n = j.height
        if n == 0:
            verdict, frac = "AMBIGUOUS_no_rsid_overlap", float("nan")
        else:
            frac = (j["POS_SUS"] == j["POS_REF"]).sum() / n
            if frac >= MATCH_38_THRESHOLD:
                verdict = "GRCh38"
            elif frac <= (1 - MATCH_38_THRESHOLD):
                verdict = "GRCh37"
            else:
                verdict = "AMBIGUOUS"
        liftover = verdict == "GRCh37"
        print(f"  overlap={n:,}  frac_pos_equal={frac:.4f}  -> {verdict}"
              f"  (liftover={'YES' if liftover else 'no'})")
        rows.append(dict(cohort=name, n_overlap=n, frac_pos_equal=round(frac, 6)
                         if n else "", verdict=verdict,
                         liftover="true" if liftover else "false"))

    out = pl.DataFrame(rows)
    out.write_csv(args.out, separator="\t")
    print(f"\nWrote {args.out}")
    print("Review verdicts, then run 02_harmonize.py (it reads this file). If any row")
    print("is AMBIGUOUS, tell me and we pick a landmark-SNP fallback before lifting.")


if __name__ == "__main__":
    sys.exit(main())
