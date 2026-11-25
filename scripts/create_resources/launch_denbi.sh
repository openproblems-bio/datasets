#!/bin/bash
# Launch a dataset workflow on the de.NBI cluster through Seqera, publishing the datasets it
# creates to LaminDB with nf-lamin.
#
#   launch_denbi.sh <workflow> <params_file> <prefix> [--large] [--labels a,b]
#
# <workflow> is a workflow namespace and name, e.g. `scrnaseq/load_cellxgene_census`, and
# <prefix> the key prefix of the artifacts (`scrnaseq` or `multimodal`). `--large` adds
# src/config/large_datasets.config for datasets that need more memory and time.
#
# Environment variables:
#   LAMIN_INSTANCE  instance to publish to (default: openproblems/data-dev)
#   ALLOW_INSTANCE  set to 1 to publish to any other instance
#   DRY_RUN         set to 1 to print the config and params instead of launching

set -euo pipefail

REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

if [[ $# -lt 3 ]]; then
  echo "usage: $0 <workflow> <params_file> <prefix> [--large] [--labels a,b]" >&2
  exit 1
fi
workflow=$1
params_file=$2
prefix=$3
shift 3

large=false
labels="$(basename "$workflow"),dataset_build"
while [[ $# -gt 0 ]]; do
  case $1 in
    --large) large=true; shift ;;
    --labels) labels=$2; shift 2 ;;
    *) echo "unknown argument: $1" >&2; exit 1 ;;
  esac
done

if [[ $prefix != scrnaseq && $prefix != multimodal ]]; then
  echo "prefix must be 'scrnaseq' or 'multimodal', got '$prefix'" >&2
  exit 1
fi

instance=${LAMIN_INSTANCE:-openproblems/data-dev}
if [[ $instance != openproblems/data-dev && ${ALLOW_INSTANCE:-0} != 1 ]]; then
  echo "refusing to publish to '$instance'; set ALLOW_INSTANCE=1 if you mean it" >&2
  exit 1
fi
publish_dir="lamin://$instance?prefix=$prefix"

tmp=$(mktemp -d)

# later files override earlier ones
config="$tmp/nextflow.config"
{
  cat common/nextflow_helpers/labels_denbi.config
  echo
  cat src/config/lamin_denbi.config
  echo
  echo "lamin { instance = '$instance' }"
  if [[ $large == true ]]; then
    echo
    cat src/config/large_datasets.config
  fi
} > "$config"

params="$tmp/params.yaml"
PUBLISH_DIR="$publish_dir" yq '.publish_dir = strenv(PUBLISH_DIR)' "$params_file" > "$params"

echo "workflow: $workflow"
echo "publish_dir: $publish_dir"
echo "config: $config"
echo "params: $params"

if [[ ${DRY_RUN:-0} == 1 ]]; then
  exit 0
fi

tw launch https://github.com/openproblems-bio/datasets.git \
  --revision build/main \
  --pull-latest \
  --main-script "target/nextflow/workflows/$workflow/main.nf" \
  --workspace 53907369739130 \
  --compute-env denbi_bibigrid_cpu \
  --params-file "$params" \
  --config "$config" \
  --labels "$labels" \
  --workspace-secrets LAMIN_BOT_API_KEY
