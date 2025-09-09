#!/bin/bash

# fetch cell cycle genes to ensure that these genes are always included
if [ ! -f /tmp/temp_g2m_genes_tirosh_mm.txt ]; then
  wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/g2m_genes_tirosh.txt -O /tmp/temp_g2m_genes_tirosh_mm.txt
fi
if [ ! -f /tmp/temp_s_genes_tirosh_mm.txt ]; then
  wget https://raw.githubusercontent.com/theislab/scib/c993ffd9ccc84ae0b1681928722ed21985fb91d1/scib/resources/s_genes_tirosh.txt -O /tmp/temp_s_genes_tirosh_mm.txt
fi

KEEP_FEATURES=`cat /tmp/temp_g2m_genes_tirosh_mm.txt /tmp/temp_s_genes_tirosh_mm.txt | paste -sd ":" -`

# example script to run the cellxgene census workflow
nextflow run openproblems-bio/datasets \
  -main-script target/nextflow/workflows/scrnaseq/process_cellxgene_census/main.nf \
  -c common/nextflow_helpers/labels_ci.config \
  -profile docker \
  --param_list scripts/new_datasets/params.yaml \
  --normalization_methods log_cp10k \
  --n_obs 600 \
  --n_vars 1500 \
  --output_dataset '$id/dataset.h5ad' \
  --output_meta '$id/dataset_metadata.yaml' \
  --output_state '$id/state.yaml' \
  --output_raw force_null \
  --output_normalized force_null \
  --output_pca force_null \
  --output_hvg force_null \
  --output_knn force_null \
  --publish_dir foobar \
  --do_subsample true \
  --keep_features "$KEEP_FEATURES"
