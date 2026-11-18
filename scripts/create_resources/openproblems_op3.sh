#!/bin/bash
# Load the OP3 dataset as a raw dataset.
set -e

params_file="/tmp/datasets_op3.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_op3
    input: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad
    dataset_name: "OP3: single-cell multimodal dataset in PBMCs for perturbation prediction benchmarking"
    dataset_summary: "The Open Problems Perurbation Prediction (OP3) dataset with small molecule perturbations in PBMCs"
    dataset_description: "The OP3 dataset is to-date the largest single-cell small molecule perturbation dataset in primary tissue with multiple donor replicates."
    dataset_url: "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE279nnn/GSE279945/suppl/GSE279945_sc_counts_processed.h5ad"
    dataset_reference: GSE279945
do_subsample: False
output_raw: '$id/raw.anndata.zarr'
output_state: '$id/state.yaml'
publish_dir: s3://openproblems-data/resources/datasets
HERE

cat > "/tmp/nextflow.config" << 'HERE'
process {
  withName:'.*publishStatesProc' {
    memory = '100GB'
    disk = '1000GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/datasets.git \
 --revision build/main \
 --pull-latest \
 --main-script target/nextflow/workflows/scrnaseq/load_openproblems_op3/main.nf \
 --workspace 53907369739130 \
 --params-file "$params_file" \
 --labels op3,dataset_loader \
 --config /tmp/nextflow.config

