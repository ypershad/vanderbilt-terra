# CLAUDE.md — Vanderbilt Terra workspace (single-cell, mutational signatures, GWAS meta)

> Read `_shared/CLAUDE-base.md` first. Rules below override the base on conflict.

## What this is

The flexible workspace. Data lives in the **Google bucket** — self-processed,
downloaded, or generated. **No BigQuery, no BioVU/AoU/biobank query layer.**
UK Biobank here is **downloaded data, not the RAP**.

Three workstreams: **single-cell analysis**, **mutational signatures**
(bcftools/samtools → SigProfiler), and **GWAS meta-analysis** (REGENIE → METAL).
Both WDL workflows and interactive notebooks.

---

## Data access

- Bucket-driven: `gs://fc-secure-827e3d9c-02b8-499f-91af-c2142b7d2074` for temporary files OR `gs://bicklab-main-storage/Users/Yash_Pershad` for more permanent files or sharing with others.
- UKB (downloaded, not RAP): <!-- FILL IN layout — which fields, formats, dirs -->
- Own sequencing: <!-- FILL IN format (CRAM/VCF/gVCF), naming, staging location -->
- Single-cell inputs: <!-- FILL IN (10x mtx? h5ad? Seurat .rds?) -->

**Path rules (learned the hard way):**
- Never hardcode a bucket path inside a WDL task — pass it as an input.
- In notebooks, anchor everything to one `DATA_DIR` constant + a `p()` helper.
  Bare relative paths break when the notebook launches from a different directory.
  ```r
  DATA_DIR <- "/home/jupyter/bicklab-pershad/edit/<subdir>"   # FILL IN subdir
  DATA_DIR <- "/home/jupyter/<workspace>/edit/<subdir>"   # FILL IN
  p <- function(...) file.path(DATA_DIR, ...)
  ```

---

## GWAS meta-analysis (METAL) — CONFIRMED conventions

**Install:** METAL from https://github.com/statgen/METAL (source:
https://csg.sph.umich.edu/abecasis/Metal/download/). Invoked as
`./generic-metal/executables/metal <script>.txt`. <!-- CONFIRM path in this workspace -->

**Standard script header (my default):**
```
SCHEME STDERR
COLUMNCOUNTING LENIENT
FORCEUPPER ON
STRICT ON
```
For multi-cohort work I also use, as needed:
`GENOMICCONTROL ON|OFF` (per-cohort), `AVERAGEFREQ ON`, `MINMAXFREQ ON`,
`CUSTOMVARIABLE N`, `SEPARATOR TAB|WHITESPACE`, `VERBOSE ON`.

**Per-cohort blocks.** METAL settings are stateful — labels persist until changed,
so each cohort gets its own block before its `PROCESS`. Two labeling styles both
appear in my scripts; keep them consistent within a script:
- Short form: `MARKER ID` / `ALLELE ALLELE0 ALLELE1` / `EFFECT BETA` / `STDERR SE` / `PVALUE PVAL`
- Long form: `MARKERLABEL` / `ALLELELABELS` / `EFFECTLABEL` / `STDERRLABEL` / `FREQLABEL`

**Filters** (applied per-cohort, then cleared with `REMOVE FILTERS`):
```
ADDFILTER eMAC >= 100
ADDFILTER info > 0.4
ADDFILTER HWE_P > 0.00000005
REMOVE FILTERS
```

**Close every script with:**
```
OUTFILE <pheno>_meta_results .tbl
ANALYZE HETEROGENEITY
QUIT
```

**Harmonization before METAL.** Cohorts arrive in different formats; normalize each
to a common schema first (see `scripts/chip_gwas_metal.sh` for the working example):
`CHR POS ID REF ALT FREQ BETA SE PVAL N`
- Standard variant ID: `chr{CHR}:{POS}:{REF}:{ALT}`
- Strip `chr` prefixes where a cohort includes them; keep ID construction consistent.
- Several sources give `LOG10P` not `PVAL` → `pval = 10^(-log10p)`.
- UKB summary files may bury effects in an INFO-style field
  (`REGENIE_BETA=`, `REGENIE_SE=`) that must be parsed out.
- Sample size is often unavailable → `N = NA`; note it rather than inventing one.
- QC after formatting: count missing BETA/SE, count variants per cohort, confirm
  column counts match the header.

**Reading results:** `<pheno>_meta_results1.tbl`. Genome-wide significance
`p < 5e-8`; when the table carries `-log10P`, that's `> 7.3`.

**Summary-stat I/O:** use **polars** in Python to read GWAS summary stats (much
faster than pandas at this size). Write the minimal columns needed downstream.

**Plots:** `qqman` in R (`manhattan()`, `qq()`) from a minimal TSV — typically just
`SNP CHR BP P`. Read with `data.table::fread`, coerce `CHR`/`BP` to numeric before
plotting. Always report lambda alongside the QQ.

**Upstream GWAS:** REGENIE <!-- FILL IN standard step1/step2 flags -->

---

## Mutational signatures

- **SigProfiler** (MatrixGenerator → Extractor / Assignment)
  <!-- FILL IN components, versions, install method (pip? conda?) -->
- Reference build: <!-- FILL IN hg19/hg38 — MUST match the VCFs; state explicitly -->
- COSMIC signature set/version: <!-- FILL IN e.g. v3.4 -->
- Pipeline: VCF → matrix generation → extraction or assignment → plots
- Install/genome-download step is slow and one-time — keep it in a separate,
  clearly-marked setup cell, never inline with analysis.

## Variant tooling

- **bcftools / samtools** <!-- FILL IN versions -->
- Standard filters/flags: <!-- FILL IN e.g. bcftools view -f PASS, -i 'INFO/AF<0.01' -->
- Always `tabix`/index after writing a VCF; use `--threads` where supported.
  <!-- FILL IN your default thread count -->

## Single-cell

<!-- FILL IN which stack is canonical -->
- R: Seurat <!-- version -->; Python: scanpy/anndata <!-- if used -->
- Integration/batch correction: <!-- Harmony? scVI? -->
- Standard workflow: QC (nFeature / nCount / percent.mt cutoffs) → normalize →
  HVG → PCA → integrate → neighbors/UMAP → cluster → markers
  <!-- FILL IN your standard cutoffs and dims so they stop being re-invented -->

---

## Compute patterns

- **WDL (Cromwell):** anything scaled across samples or chromosomes.
  Default runtime: <!-- FILL IN docker image, memory, cpu, disks -->
  One task per logical step, explicit `runtime{}`, all paths as inputs.
- **Notebooks:** interactive analysis and plotting downstream of WDL/CLI outputs.
- Long-running CLI jobs (METAL, SigProfiler, bcftools over many files): wrap in a
  shell script under `scripts/`, run with `nohup ... &`, check with `jobs -l`.

## Terra launch constraint (same as our BioVU workspace)

- Terra's Jupyter GUI only opens notebooks located directly in the `edit/` root.
- Repo lives under `edit/`: ~/bicklab-pershad/edit/vanderbilt-terra
  - `./bin/nb-open.sh <name>` — materialize `edit/<name>.ipynb` from `notebooks/<name>.R|.py`
  - `./bin/nb-save.sh <name> ["msg"]` — sync edits back to the repo script, commit, push
- Edit the `.R`/`.py`, never the `.ipynb`. Outputs stripped; `.ipynb` gitignored.

## Style

- Concise, correct, minimal. Follow instructions; don't over-engineer.
- Complex task → reason step by step first, then write code.
- Comment the *why* (why this filter, why this scheme), not the obvious *what*.
- Don't import libraries not in the image without flagging it.
- Keep one constants block per notebook (DATA_DIR, bucket, build) — no scattered literals.
