#!/usr/bin/env python3
"""
00_inspect_cohorts.py — DISCOVERY ONLY. Hardcodes NO column mappings.

Purpose: for each of the 13 hematuria GWAS summary-stat files, print the header
and first 3 data rows so column meanings can be confirmed BEFORE any parser is
written (Task 1, step 1). Also emits cheap format/build *hints* (never a decision).

Why this exists as its own step: per _shared/CLAUDE-base.md I write code on the
laptop and run it on the Terra VM. This script is the thing you run first on the
VM; paste its output back and the per-cohort parser configs get written with
certainty instead of guesses.

Design notes:
  - Streams only the first N lines of each gzip (no full-file load — you never need
    to decompress 3.8 GB to read a header). This is the ONE place we don't follow
    the "load everything into RAM" constraint, because it's pure discovery.
  - Prints summary stats only (aggregate columns + a few rows). No individual-level
    data, no identifiers — safe under the controlled-access rules.
  - Delimiter is sniffed, not assumed. Comment/meta lines (e.g. FinnGen '#chrom')
    are surfaced, not silently skipped.

Usage on the VM:
    python3 scripts/00_inspect_cohorts.py \
        --data-dir hematuria_gwas/biovu_results \
        --n-rows 3 | tee hematuria_gwas/inspect_report.txt

Then paste inspect_report.txt back into the chat.
"""

import argparse
import gzip
import os
import sys
from collections import Counter

# The 13 expected files. Metadata here is only what the FILENAMES / your brief
# assert — NOT column knowledge. `build_claim` is what the source says; it still
# gets confirmed from the header/rows below, and Taiwan is explicitly unknown.
EXPECTED = [
    # BioVU AGD250k, REGENIE 4.0 output. GRCh38 per brief. 2 ancestries x 4 phenodefs.
    ("agd250k_eur.593.pvalue_ordered.regenie.gz",      dict(cohort="biovu_agd250k", ancestry="EUR", phenodef="593",      build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_eur.593.2.pvalue_ordered.regenie.gz",    dict(cohort="biovu_agd250k", ancestry="EUR", phenodef="593.2",    build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_eur.GU_593.pvalue_ordered.regenie.gz",   dict(cohort="biovu_agd250k", ancestry="EUR", phenodef="GU_593",   build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_eur.GU_593.2.pvalue_ordered.regenie.gz", dict(cohort="biovu_agd250k", ancestry="EUR", phenodef="GU_593.2", build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_afr.593.pvalue_ordered.regenie.gz",      dict(cohort="biovu_agd250k", ancestry="AFR", phenodef="593",      build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_afr.593.2.pvalue_ordered.regenie.gz",    dict(cohort="biovu_agd250k", ancestry="AFR", phenodef="593.2",    build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_afr.GU_593.pvalue_ordered.regenie.gz",   dict(cohort="biovu_agd250k", ancestry="AFR", phenodef="GU_593",   build_claim="GRCh38", fmt_hint="regenie")),
    ("agd250k_afr.GU_593.2.pvalue_ordered.regenie.gz", dict(cohort="biovu_agd250k", ancestry="AFR", phenodef="GU_593.2", build_claim="GRCh38", fmt_hint="regenie")),
    # External cohorts
    ("mgi_eur_X593.output_GRCh38.gz",                  dict(cohort="mgi",     ancestry="EUR", phenodef="X593",   build_claim="GRCh38 (filename)", fmt_hint="unknown")),
    ("ukb_pheweb_phenocode-593.tsv.gz",                dict(cohort="ukb",     ancestry="EUR", phenodef="593",    build_claim="GRCh37? (PheWeb — CONFIRM)", fmt_hint="pheweb?")),
    ("summary_stats_finngen_R7_R18_UNSPE_HAEMATU.gz",  dict(cohort="finngen", ancestry="FIN", phenodef="R18_UNSPE_HAEMATU", build_claim="GRCh38 (R7)", fmt_hint="finngen")),
    ("593_Taiwan_biobank_Aug132025.tsv.gz",            dict(cohort="taiwan",  ancestry="EAS", phenodef="593",    build_claim="UNKNOWN — CHECK EXPLICITLY", fmt_hint="unknown")),
    ("phenocode-KoGES_BLOODU.tsv.gz",                  dict(cohort="koges",   ancestry="EAS", phenodef="BLOODU", build_claim="GRCh37? (CONFIRM)", fmt_hint="pheweb?")),
]

# Header tokens that *hint* at a format. Hints only — the human confirms.
FORMAT_SIGNATURES = {
    "regenie": {"CHROM", "GENPOS", "ALLELE0", "ALLELE1", "A1FREQ", "LOG10P"},
    "finngen/pheweb": {"#chrom", "mlogp", "af_alt", "nearest_genes", "rsids"},
    "saige":   {"AC_Allele2", "AF_Allele2", "Tstat", "var", "p.value"},
    "plink2":  {"OBS_CT", "A1_FREQ", "LOG10_P", "OR"},
}

# Column-name tokens that hint at genome build. NEVER decisive — a real build
# check compares a few known-position variants; flagged for human follow-up.
BUILD_HINT_TOKENS = {"build", "genome", "grch", "hg19", "hg38", "assembly"}


def sniff_delimiter(line):
    """Pick the delimiter that splits into the most fields. Returns (name, char)."""
    candidates = [("tab", "\t"), ("comma", ","), ("space", " ")]
    best = max(candidates, key=lambda c: len(line.split(c[1])))
    # whitespace-collapsed fallback (REGENIE/METAL often space-padded)
    ws_fields = len(line.split())
    if ws_fields > len(line.split(best[1])):
        return ("whitespace", None)
    return best


def split_line(line, delim_char):
    return line.split() if delim_char is None else line.split(delim_char)


def read_head(path, n_rows, max_scan=50):
    """Return (raw_lines, comment_lines) reading at most a few lines. Streams gzip."""
    raw, comments = [], []
    with gzip.open(path, "rt", errors="replace") as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if i >= max_scan:
                break
            # capture leading pure-comment/meta lines that are NOT the header
            if line.startswith("##") or (line.startswith("#") and not raw and _looks_like_meta(line)):
                comments.append(line)
                continue
            raw.append(line)
            if len(raw) >= n_rows + 1:  # header + n_rows
                break
    return raw, comments


def _looks_like_meta(line):
    # '#chrom\t...' is a header (many tab fields); '## key=value' or '# generated by' is meta.
    return len(line.split()) <= 3 or "=" in line


def format_hint(header_tokens):
    hits = []
    toks = set(header_tokens) | {t.lower() for t in header_tokens}
    for fmt, sig in FORMAT_SIGNATURES.items():
        sig_l = {s.lower() for s in sig}
        overlap = len(sig_l & {t.lower() for t in header_tokens})
        if overlap >= 2:
            hits.append(f"{fmt}(+{overlap})")
    return ", ".join(hits) if hits else "no strong signature"


def build_hint(header_tokens):
    flagged = [t for t in header_tokens if any(k in t.lower() for k in BUILD_HINT_TOKENS)]
    return flagged


def inspect_one(path, meta, n_rows):
    print("=" * 88)
    fname = os.path.basename(path)
    print(f"FILE: {fname}")
    print(f"  cohort={meta['cohort']}  ancestry={meta['ancestry']}  "
          f"phenodef={meta['phenodef']}")
    print(f"  build (claimed): {meta['build_claim']}   fmt (guess): {meta['fmt_hint']}")

    if not os.path.exists(path):
        print("  !! FILE NOT FOUND at this path — check --data-dir")
        print()
        return

    size_gb = os.path.getsize(path) / 1e9
    print(f"  size: {size_gb:.2f} GB gzipped")
    # note tabix sibling (position-sorted / bgzipped) if present
    if os.path.exists(path + ".tbi"):
        print("  note: .tbi present -> bgzipped + tabix-indexed (position-sorted)")

    try:
        raw, comments = read_head(path, n_rows)
    except Exception as e:
        print(f"  !! failed to read: {e!r}")
        print()
        return

    if comments:
        print(f"  leading meta/comment lines ({len(comments)}):")
        for c in comments[:5]:
            print(f"    | {c[:120]}")

    if not raw:
        print("  !! no data lines found")
        print()
        return

    header = raw[0]
    dname, dchar = sniff_delimiter(header)
    cols = split_line(header, dchar)
    print(f"  delimiter (sniffed): {dname}")
    print(f"  n_columns: {len(cols)}")
    print(f"  raw header line:\n    {header[:400]}")
    print("  columns (index: name):")
    for i, c in enumerate(cols):
        print(f"    [{i:>2}] {c}")

    print(f"  format signature match: {format_hint(cols)}")
    bh = build_hint(cols)
    if bh:
        print(f"  !! build-related column(s) present -> read their VALUE below: {bh}")

    print(f"  first {n_rows} data rows (as index:name=value):")
    for r_i, row in enumerate(raw[1:1 + n_rows], start=1):
        vals = split_line(row, dchar)
        pairs = [f"{cols[i] if i < len(cols) else f'col{i}'}={v}" for i, v in enumerate(vals)]
        # show first 14 fields inline; note if truncated
        shown = "  ".join(pairs[:14])
        more = f"  ... (+{len(vals) - 14} more fields)" if len(vals) > 14 else ""
        if len(vals) != len(cols):
            more += f"  [!! {len(vals)} values vs {len(cols)} header cols]"
        print(f"    row{r_i}: {shown}{more}")
    print()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", default="../hematuria_gwas/biovu_results",
                    help="dir holding the 13 .gz files. Default assumes you run from "
                         "the repo (edit/vanderbilt-terra); data is the sibling "
                         "edit/hematuria_gwas/biovu_results.")
    ap.add_argument("--n-rows", type=int, default=3, help="data rows to print per file")
    args = ap.parse_args()

    print(f"Inspecting {len(EXPECTED)} expected files under: {args.data_dir}\n")

    # First: reconcile expected vs actual directory contents so nothing is missed.
    if os.path.isdir(args.data_dir):
        actual = {f for f in os.listdir(args.data_dir)
                  if f.endswith((".gz", ".tsv", ".txt", ".regenie"))}
        expected_names = {f for f, _ in EXPECTED}
        missing = expected_names - actual
        extra = actual - expected_names - {n + ".tbi" for n in expected_names}
        if missing:
            print(f"!! MISSING expected files: {sorted(missing)}\n")
        if extra:
            print(f"?? EXTRA files in dir (not in expected list): {sorted(extra)}\n")
    else:
        print(f"!! --data-dir does not exist or is not a directory: {args.data_dir}\n")

    for fname, meta in EXPECTED:
        inspect_one(os.path.join(args.data_dir, fname), meta, args.n_rows)

    print("=" * 88)
    print("DONE. Paste this whole report back so parser configs can be written with")
    print("certainty. Note especially: Taiwan build, UKB/KoGES build, any OR-vs-BETA")
    print("columns, and whether REGENIE files carry a Firth TEST column.")


if __name__ == "__main__":
    sys.exit(main())
