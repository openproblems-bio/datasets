#!/bin/bash
# Tests for launch_denbi.sh, run without launching anything (DRY_RUN=1).
#
#   bash scripts/create_resources/test_launch_denbi.sh

set -euo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

launcher=scripts/create_resources/launch_denbi.sh
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

cat > "$tmp/params.yaml" << 'HERE'
param_list:
  - id: foo
publish_dir: s3://somewhere/else
HERE

failures=0
fail() { echo "FAIL: $*"; failures=$((failures + 1)); }
pass() { echo "ok: $*"; }

# extract the value of a 'key: value' line from the dry-run output
field() { sed -n "s/^$1: //p" <<< "$2"; }

# the merged config contains labels_denbi.config, then lamin_denbi.config, then the instance
out=$(DRY_RUN=1 bash "$launcher" scrnaseq/load_cellxgene_census "$tmp/params.yaml" scrnaseq)
config=$(field config "$out")
params=$(field params "$out")
line_of() { grep -n -F -- "$1" "$config" | head -1 | cut -d: -f1; }
l_denbi=$(line_of "executor = 'slurm'")
l_lamin=$(line_of "nf-lamin@")
l_instance=$(line_of "instance = 'openproblems/data-dev'")
if [[ -n $l_denbi && -n $l_lamin && -n $l_instance && $l_denbi -lt $l_lamin && $l_lamin -lt $l_instance ]]; then
  pass "merged config order"
else
  fail "merged config order (denbi=$l_denbi lamin=$l_lamin instance=$l_instance)"
fi
if grep -q "450.GB" "$config"; then fail "large overlay without --large"; else pass "no large overlay by default"; fi

# the params point at the Lamin instance, and the original publish_dir is gone
if [[ $(yq '.publish_dir' "$params") == "lamin://openproblems/data-dev?prefix=scrnaseq" ]]; then
  pass "publish_dir"
else
  fail "publish_dir is '$(yq '.publish_dir' "$params")'"
fi
if [[ $(yq '.param_list[0].id' "$params") == "foo" ]]; then pass "params kept"; else fail "params lost"; fi

# --large appends large_datasets.config, after the instance
out=$(DRY_RUN=1 bash "$launcher" scrnaseq/load_cellxgene_census "$tmp/params.yaml" scrnaseq --large)
config=$(field config "$out")
if grep -q "450.GB" "$config"; then pass "--large overlay"; else fail "--large overlay missing"; fi

# any other instance is refused, unless explicitly allowed
if LAMIN_INSTANCE=openproblems/data DRY_RUN=1 bash "$launcher" scrnaseq/load_cellxgene_census "$tmp/params.yaml" scrnaseq > /dev/null 2>&1; then
  fail "openproblems/data accepted without ALLOW_INSTANCE"
else
  pass "openproblems/data refused"
fi
out=$(ALLOW_INSTANCE=1 LAMIN_INSTANCE=openproblems/data DRY_RUN=1 bash "$launcher" scrnaseq/load_cellxgene_census "$tmp/params.yaml" multimodal)
if [[ $(yq '.publish_dir' "$(field params "$out")") == "lamin://openproblems/data?prefix=multimodal" ]]; then
  pass "ALLOW_INSTANCE"
else
  fail "ALLOW_INSTANCE"
fi

# an unknown prefix is refused
if DRY_RUN=1 bash "$launcher" scrnaseq/load_cellxgene_census "$tmp/params.yaml" spatial > /dev/null 2>&1; then
  fail "unknown prefix accepted"
else
  pass "unknown prefix refused"
fi

[[ $failures -eq 0 ]] && echo "all tests passed" || { echo "$failures test(s) failed"; exit 1; }
