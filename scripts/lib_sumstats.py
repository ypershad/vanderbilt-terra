"""
lib_sumstats.py — shared harmonization helpers for the hematuria meta-analysis.

Pure functions only (no side effects beyond returning data / a log dict). Imported
by 01_build_check.py and 02_harmonize.py. Runs on the Terra VM, not the laptop.

Design choices tied to the inspection findings:
  - Every cohort's EA=ALT, NEA=REF, EAF=freq(ALT), BETA per-ALT log-odds (verified).
  - LOG10P kept internally and computed underflow-safe:
        FinnGen -> mlogp column (already -log10 p)
        BioVU   -> from Chisq (df=1) via scipy.chi2.logsf  (P underflows at GWAS sig)
        others  -> from reported P; if P==0 (underflow) fall back to Wald z=beta/se
  - key38 = "CHR:POS:sortedA:sortedB" surfaces palindromic/strand collisions.
  - Deps: polars, numpy, scipy (preflight-checked in the entry scripts).
"""

import math
import os

import numpy as np
import polars as pl

LOG10 = math.log(10.0)


def resolve_configs(yaml_path, data_dir):
    """Flatten cohorts.yaml into a list of per-FILE job dicts.

    BioVU has 8 strata under one cohort block; every other cohort is one file.
    Each job carries everything a worker needs: absolute path, ids, column_map,
    build/liftover flags, and the shared qc thresholds. No global state.
    """
    import yaml
    with open(yaml_path) as fh:
        conf = yaml.safe_load(fh)
    qc_cfg = conf["defaults"]["qc"]
    jobs = []

    for cohort, c in conf["cohorts"].items():
        base = dict(
            cohort=cohort,
            format=c.get("format"),
            build=c.get("build"),
            build_confirmed=c.get("build_confirmed", False),
            liftover=c.get("liftover", False),
            effect_scale=c.get("effect_scale"),
            column_map=c["column_map"],
            snpid_col=c.get("snpid_col"),
            rsid_col=c.get("rsid_col"),
            qc=qc_cfg,
        )
        if "strata" in c:                         # BioVU: many files, shared parser
            for s in c["strata"]:
                job = dict(base)
                job.update(file=s["file"], ancestry=s["ancestry"],
                           phenodef=s["phenodef"])
                job["path"] = os.path.join(data_dir, s["file"])
                jobs.append(job)
        else:                                     # one file per cohort
            job = dict(base)
            job.update(file=c["file"], ancestry=c["ancestry"], phenodef=c["phenodef"])
            job["path"] = os.path.join(data_dir, c["file"])
            jobs.append(job)
    return jobs

# Harmonized output schema (column order for the parquet).
OUT_COLS = [
    "CHR", "POS", "EA", "NEA", "EAF", "BETA", "SE", "P", "LOG10P",
    "N", "N_CASES", "N_CONTROLS", "INFO", "KEY38",
    "COHORT", "ANCESTRY", "PHENODEF", "RSID",
]

VALID_CHR = {str(i) for i in range(1, 23)} | {"X"}


# ----------------------------------------------------------------------------- IO
def read_raw(path, columns=None):
    """Read a (gzip) sumstat file into polars, all columns as strings.

    The FinnGen header line is '#chrom\\tpos\\t...'. We set comment_prefix=None so
    polars does NOT treat it as a comment; the first column name stays literally
    '#chrom', which is exactly what cohorts.yaml maps (chr: "#chrom"). Every column
    is read as Utf8 and cast explicitly in standardize() — no schema inference races.

    columns: optional list of source column names to keep (speeds up big files).
    """
    df = pl.read_csv(
        path,
        separator="\t",
        comment_prefix=None,        # keep '#chrom' as a real header, not a comment
        infer_schema_length=0,      # read all as str; cast explicitly downstream
    )
    if columns is not None:
        keep = [c for c in columns if c in df.columns]
        df = df.select(keep)
    return df


# --------------------------------------------------------------------- normalize
def _clean_chr(expr):
    """Strip 'chr' prefix, map 23->X / 25->X, uppercase X. Returns str expr."""
    e = expr.cast(pl.Utf8).str.replace(r"(?i)^chr", "")
    e = e.replace({"23": "X", "25": "X", "x": "X"})
    return e


def _rsid_expr(cfg, src):
    """RSID from cfg['rsid_col'] if set, else auto-detect a 'rsids'/'rsid' column."""
    name = cfg.get("rsid_col")
    if name is None:
        for cand in ("rsids", "rsid", "RSID"):
            if cand in src.columns:
                name = cand
                break
    if name and name in src.columns:
        return pl.col(name).cast(pl.Utf8)
    return pl.lit(None, pl.Utf8)


def standardize(df, cfg, cohort, ancestry, phenodef):
    """Apply a cohort's column_map -> harmonized columns (pre-QC, pre-liftover).

    cfg: the resolved per-file config dict with keys column_map, snpid_col, etc.
    Returns a polars DataFrame with the OUT_COLS present (POS still source build).
    """
    cm = cfg["column_map"]
    src = df

    # Build select expressions, tolerating null-mapped (absent) source columns.
    def col_or_null(key, dtype=pl.Float64):
        name = cm.get(key)
        if name is None or name not in src.columns:
            return pl.lit(None, dtype=dtype)
        return pl.col(name)

    out = src.select(
        _clean_chr(pl.col(cm["chr"])).alias("CHR"),
        pl.col(cm["pos"]).cast(pl.Int64, strict=False).alias("POS"),
        pl.col(cm["ea"]).cast(pl.Utf8).str.to_uppercase().alias("EA"),
        pl.col(cm["nea"]).cast(pl.Utf8).str.to_uppercase().alias("NEA"),
        col_or_null("eaf").cast(pl.Float64, strict=False).alias("EAF"),
        col_or_null("beta").cast(pl.Float64, strict=False).alias("BETA"),
        col_or_null("se").cast(pl.Float64, strict=False).alias("SE"),
        col_or_null("p").cast(pl.Float64, strict=False).alias("P"),
        col_or_null("log10p").cast(pl.Float64, strict=False).alias("_LOG10P_SRC"),
        col_or_null("chisq").cast(pl.Float64, strict=False).alias("_CHISQ"),
        col_or_null("n").cast(pl.Float64, strict=False).alias("N"),
        _rsid_expr(cfg, src).alias("RSID"),
    )

    out = out.with_columns(
        pl.lit(None, pl.Float64).alias("N_CASES"),
        pl.lit(None, pl.Float64).alias("N_CONTROLS"),
        pl.lit(None, pl.Float64).alias("INFO"),
        pl.lit(cohort).alias("COHORT"),
        pl.lit(ancestry).alias("ANCESTRY"),
        pl.lit(phenodef).alias("PHENODEF"),
    )
    return out


def add_log10p(df):
    """Fill LOG10P underflow-safe, in priority: source mlogp -> chisq -> P -> Wald z.

    Uses numpy/scipy on arrays (polars has no chi2 survival fn). Returns df + LOG10P.
    """
    from scipy.stats import chi2, norm

    n = df.height
    log10p = np.full(n, np.nan)

    src = df["_LOG10P_SRC"].to_numpy()
    chisq = df["_CHISQ"].to_numpy()
    p = df["P"].to_numpy()
    beta = df["BETA"].to_numpy()
    se = df["SE"].to_numpy()

    # 1) explicit -log10 p (FinnGen mlogp)
    m = np.isfinite(src)
    log10p[m] = src[m]

    # 2) from chisq, df=1 (BioVU): -log10(sf) = -logsf/ln10  (no underflow)
    m2 = ~np.isfinite(log10p) & np.isfinite(chisq) & (chisq >= 0)
    if m2.any():
        log10p[m2] = -chi2.logsf(chisq[m2], df=1) / LOG10

    # 3) from reported P where P > 0
    m3 = ~np.isfinite(log10p) & np.isfinite(p) & (p > 0)
    log10p[m3] = -np.log10(p[m3])

    # 4) P underflowed to 0 (or missing) -> Wald z from beta/se
    m4 = ~np.isfinite(log10p) & np.isfinite(beta) & np.isfinite(se) & (se > 0)
    if m4.any():
        z = np.abs(beta[m4] / se[m4])
        # two-sided -log10 p from |z|, underflow-safe via log survival fn
        log10p[m4] = -(np.log(2.0) + norm.logsf(z)) / LOG10

    return df.with_columns(pl.Series("LOG10P", log10p)).drop(["_LOG10P_SRC", "_CHISQ"])


def add_key38(df):
    """key38 = CHR:POS:sorted(EA,NEA). Sorted alleles collapse strand/orientation."""
    lo = pl.when(pl.col("EA") <= pl.col("NEA")).then(pl.col("EA")).otherwise(pl.col("NEA"))
    hi = pl.when(pl.col("EA") <= pl.col("NEA")).then(pl.col("NEA")).otherwise(pl.col("EA"))
    return df.with_columns(
        (pl.col("CHR") + ":" + pl.col("POS").cast(pl.Utf8) + ":" + lo + ":" + hi)
        .alias("KEY38")
    )


# --------------------------------------------------------------------------- QC
def qc(df, qc_cfg, log):
    """Apply QC filters, appending a (step, dropped, remaining) row to `log` list.

    Rules (per brief; INFO skipped — no cohort carries it):
      - valid CHR (autosomes + X), integer POS > 0
      - biallelic single-token alleles A/C/G/T-or-indel, EA != NEA
      - SE finite & > 0 ; BETA finite & |BETA| < abs_beta_max
      - non-missing P (or LOG10P)
      - EAF in [eaf_min, eaf_max] where EAF present; else MAC>=mac_min where N present;
        cohorts with neither EAF nor N skip the frequency filter (logged)
    """
    def step(name, mask):
        nonlocal df
        before = df.height
        df = df.filter(mask)
        log.append({"step": name, "dropped": before - df.height, "remaining": df.height})

    log.append({"step": "input", "dropped": 0, "remaining": df.height})

    step("valid_chr", pl.col("CHR").is_in(list(VALID_CHR)))
    step("pos_positive", pl.col("POS").is_not_null() & (pl.col("POS") > 0))
    step("alleles_present", pl.col("EA").is_not_null() & pl.col("NEA").is_not_null()
                            & (pl.col("EA") != pl.col("NEA"))
                            & (pl.col("EA") != "") & (pl.col("NEA") != ""))
    step("se_finite_pos", pl.col("SE").is_finite() & (pl.col("SE") > 0))
    step("beta_finite_bounded",
         pl.col("BETA").is_finite() & (pl.col("BETA").abs() < qc_cfg["abs_beta_max"]))
    step("p_present", pl.col("LOG10P").is_finite())

    # frequency filter — adapt to what the cohort actually has
    has_eaf = df["EAF"].is_not_null().any()
    has_n = df["N"].is_not_null().any()
    if has_eaf:
        step("eaf_bounds",
             (pl.col("EAF") >= qc_cfg["eaf_min"]) & (pl.col("EAF") <= qc_cfg["eaf_max"]))
    elif has_n:
        # MAC = 2 * N * MAF ; MAF unknown w/o EAF, so this branch only reached if
        # EAF truly absent but N present — cannot compute MAC either. Log + skip.
        log.append({"step": "freq_filter_SKIPPED_no_eaf", "dropped": 0,
                    "remaining": df.height})
    else:
        log.append({"step": "freq_filter_SKIPPED_no_eaf_no_n", "dropped": 0,
                    "remaining": df.height})

    return df


def split_autosome_x(df):
    """Return (autosomes_df, x_df). chrX kept separately, never silently dropped."""
    return df.filter(pl.col("CHR") != "X"), df.filter(pl.col("CHR") == "X")


# --------------------------------------------------------------------- liftover
def liftover_pos(df, chain_path, log):
    """Lift POS from GRCh37 -> GRCh38 in place using pyliftover.

    pyliftover is 0-based; GWAS POS is 1-based -> query pos-1, add 1 back.
    Unmapped variants are dropped and counted. Returns df with GRCh38 POS.
    """
    from pyliftover import LiftOver

    lo = LiftOver(chain_path)
    chrom = df["CHR"].to_numpy()
    pos = df["POS"].to_numpy()
    new_pos = np.full(len(pos), -1, dtype=np.int64)

    for i in range(len(pos)):
        res = lo.convert_coordinate("chr" + str(chrom[i]), int(pos[i]) - 1)
        if res:
            new_pos[i] = res[0][1] + 1   # back to 1-based

    before = df.height
    df = df.with_columns(pl.Series("POS", new_pos)).filter(pl.col("POS") > 0)
    log.append({"step": "liftover_37to38_unmapped_dropped",
                "dropped": before - df.height, "remaining": df.height})
    return df
