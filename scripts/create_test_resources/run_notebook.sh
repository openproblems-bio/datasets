#!/bin/bash

# Run a test resource notebook and save it, including its report, to LaminDB.
#
# Usage:
#   uv run --with 'lamindb[aws]' --with nbconvert --with ipykernel \
#     scripts/create_test_resources/run_notebook.sh scripts/create_test_resources/pancreas.ipynb

set -eo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
notebook=$(realpath --relative-to="$REPO_ROOT" "$1")
cd "$REPO_ROOT"

# the .lamin dev-dir in the repo root selects the instance
instance=$(python -c "import lamindb_setup as ls; print(ls.settings.instance.slug)")
if [ "$instance" != "openproblems/data-dev" ]; then
  echo "Connected to '$instance' instead of 'openproblems/data-dev'" >&2
  exit 1
fi

# nbconvert writes the outputs into the notebook, restore the committed version afterwards
backup=$(mktemp)
cp "$notebook" "$backup"
trap 'cp "$backup" "$notebook"; rm -f "$backup"' EXIT

jupyter nbconvert --to notebook --execute --inplace "$notebook" --ExecutePreprocessor.timeout=-1

# ln.finish() doesn't save a report for nbconvert runs, lamin save does
lamin save "$notebook"
