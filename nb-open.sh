#!/usr/bin/env bash
# nb-open.sh — materialize a repo notebook script as a runnable .ipynb in Terra's edit/ root.
# Terra's Jupyter GUI only launches notebooks from edit/, so the .ipynb must live
# there even though its source (.R/.py) lives in the repo.
#
# Usage:  ./bin/nb-open.sh <notebook_name_without_extension>
#
# CONFIG: set EDIT_DIR once for this workspace (see FILL IN below), or export
# TERRA_EDIT_DIR in your shell to override without editing this file.

set -euo pipefail

# <<< FILL IN: this workspace's Terra edit/ directory >>>
EDIT_DIR="${TERRA_EDIT_DIR:-$HOME/bicklab-pershad/edit}"
REPO_DIR="${TERRA_REPO_DIR:-$EDIT_DIR/vanderbilt-terra}"

NAME="${1:?Usage: nb-open.sh <notebook_name_without_extension>}"

# Find the source script: .R or .py
SRC=""
for ext in R py; do
  [[ -f "$REPO_DIR/notebooks/$NAME.$ext" ]] && SRC="$REPO_DIR/notebooks/$NAME.$ext" && break
done
[[ -n "$SRC" ]] || { echo "ERROR: no notebooks/$NAME.R or .py in $REPO_DIR" >&2; exit 1; }

DEST="$EDIT_DIR/$NAME.ipynb"

if [[ -f "$DEST" ]]; then
  echo "WARNING: $DEST already exists."
  echo "If it has unsaved edits, run nb-save.sh first, or delete it to regenerate."
  read -r -p "Overwrite $DEST? [y/N] " ans
  [[ "$ans" == "y" || "$ans" == "Y" ]] || { echo "Aborted."; exit 1; }
fi

jupytext --to notebook "$SRC" -o "$DEST"
echo "Opened notebook at: $DEST"
echo "-> open '$NAME.ipynb' from the Jupyter file browser (it's in edit/)."
