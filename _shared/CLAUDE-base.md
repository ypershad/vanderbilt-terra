# CLAUDE-base.md — shared conventions (all workspaces)

This file holds rules common to every Terra/AoU workspace repo. Each repo's
`CLAUDE.md` says "read this first," then adds workspace-specific rules that
**override** anything here on conflict.

## Execution model

- I write/review code on my laptop with Claude Code; **I run everything myself**
  on the cloud workspace. Do not assume local execution, local filesystem, or a
  laptop GPU. No code should auto-run or call the workspace API/compute.
- Code reaches the workspace via `git pull`; live tweaks come back via `git push`.

## Hard rules — controlled-access data

- All of these workspaces touch controlled-access data. **Never** write code that
  egresses, downloads to a personal machine, emails, or otherwise moves
  individual-level data outside the secure workspace.
- Never hardcode credentials, tokens, service-account keys, or auth output. Auth is
  handled by the workspace runtime.
- Never emit patient/sample identifiers into committed logs, notebook output,
  commit messages, or filenames.
- Reference bucket paths / table names freely; row-level data, never.

## Notebooks

- **The source of truth in git is a text script, never the `.ipynb`.** The `.ipynb`
  is *generated on Terra* from the script via Jupytext, and is gitignored.
- Pairing is percent format, per language:
  - **R** notebooks ↔ **`.R`** (`# %% [markdown]` / `# %%` cell markers)
  - **Python** notebooks ↔ **`.py`** (same markers)
- **Claude edits the `.R` / `.py`, never the `.ipynb` JSON.**
- Terra workflow: `git pull` → `jupytext --to notebook <file>.R` to materialize the
  `.ipynb` → run it. After a live edit on Terra, `jupytext --sync <file>.R` writes
  changes back into the script; commit and push the script.
- Because outputs live only in the Terra-side `.ipynb` (gitignored), no cell output,
  data, or identifiers ever land in git. GitHub shows code; Terra shows the run.

## Style

- Concise, correct, minimal. Follow my instructions; don't over-engineer.
- Complex task → reason step by step first, then write code.
- Comment the *why*, not the obvious *what*.
- Don't import libraries that aren't in that workspace's image without flagging it.

## Git hygiene

- Data files, credentials, and `.ipynb_checkpoints` are gitignored — keep it that way.
- Small, descriptive commits. No data or identifiers in commit messages.
