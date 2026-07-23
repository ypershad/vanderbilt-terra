#!/usr/bin/env bash
# nb-save.sh — sync a Terra-edited .ipynb (in edit/) back to its repo script,
# then stage/commit/push. The script (no outputs) stays git's source of truth.
#
# Usage:  ./bin/nb-save.sh <notebook_name> ["commit message"]

set -euo pipefail

# <<< FILL IN: this workspace's Terra edit/ directory >>>
EDIT_DIR="${TERRA_EDIT_DIR:-$HOME/bicklab-pershad/edit}"
REPO_DIR="${TERRA_REPO_DIR:-$EDIT_DIR/vanderbilt-terra}"

NAME="${1:?Usage: nb-save.sh <notebook_name> [commit message]}"
MSG="${2:-edits from Terra: $NAME}"

SRC="$EDIT_DIR/$NAME.ipynb"
[[ -f "$SRC" ]] || { echo "ERROR: notebook not found in edit/: $SRC" >&2; exit 1; }

# Determine target extension from the existing repo script
DEST=""; FMT=""
if [[ -f "$REPO_DIR/notebooks/$NAME.R" ]];  then DEST="$REPO_DIR/notebooks/$NAME.R";  FMT="R:percent";  fi
if [[ -f "$REPO_DIR/notebooks/$NAME.py" ]]; then DEST="$REPO_DIR/notebooks/$NAME.py"; FMT="py:percent"; fi
if [[ -z "$DEST" ]]; then
  # New notebook: infer from kernel language
  LANG=$(python3 -c "import json,sys; print(json.load(open('$SRC'))['metadata'].get('kernelspec',{}).get('language','python'))")
  if [[ "$LANG" == "R" ]]; then DEST="$REPO_DIR/notebooks/$NAME.R"; FMT="R:percent";
  else DEST="$REPO_DIR/notebooks/$NAME.py"; FMT="py:percent"; fi
fi

jupytext --to "$FMT" "$SRC" -o "$DEST"

cd "$REPO_DIR"
git add "$DEST"
if git diff --cached --quiet; then
  echo "No changes to commit for $(basename "$DEST")."
  exit 0
fi
git commit -m "$MSG"
git push
echo "Pushed $(basename "$DEST") to origin."
