#!/bin/bash

set -e

params_file="/tmp/datasets_openproblems_neurips2022_params.yaml"

cat > "$params_file" << 'HERE'
param_list:
  - id: openproblems_neurips2022/pbmc_cite
    input_mod1: s3://openproblems-nextflow/datasets_private/neurips2022/cite_rna_merged.h5ad
    input_mod2: s3://openproblems-nextflow/datasets_private/neurips2022/cite_prot_merged.h5ad
    mod1: GEX
    mod2: ADT
    dataset_name: OpenProblems NeurIPS2022 CITE-Seq
    dataset_summary: Single-cell CITE-Seq (GEX+ADT) data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells from 4 healthy human donors, over 4 days.
    dataset_description: "Single-cell CITE-Seq data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells (HSPCs) from 4 healthy human donors, measured on days 2, 3, 4 and 7 of a 10-day time course. Gene expression was measured together with the abundance of 140 cell surface proteins. The dataset was generated for the Multimodal Single-Cell Integration competition at NeurIPS 2022, in which the abundance of the proteins had to be predicted from gene expression on a later day. The cells were annotated with seven progenitor cell types."

  - id: openproblems_neurips2022/pbmc_multiome
    input_mod1: s3://openproblems-nextflow/datasets_private/neurips2022/multiome_rna_merged.h5ad
    input_mod2: s3://openproblems-nextflow/datasets_private/neurips2022/multiome_atac_merged.h5ad
    mod1: GEX
    mod2: ATAC
    dataset_name: OpenProblems NeurIPS2022 Multiome
    dataset_summary: Single-cell Multiome (GEX+ATAC) data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells from 4 healthy human donors, over 5 days.
    dataset_description: "Single-cell Multiome data of mobilised peripheral blood CD34+ hematopoietic stem and progenitor cells (HSPCs) from 4 healthy human donors, measured on days 2, 3, 4, 7 and 10 of a 10-day time course with the 10X Multiome Gene Expression and Chromatin Accessibility kit. The dataset was generated for the Multimodal Single-Cell Integration competition at NeurIPS 2022, in which gene expression had to be predicted from chromatin accessibility on a later day. The cells were annotated with seven progenitor cell types; the cell types of the cells in the private test set were not released."

dataset_url: "https://www.kaggle.com/competitions/open-problems-multimodal/data"
dataset_reference: lance2024predicting
output_raw: '$id/raw.mudata.zarr'
output_state: '$id/state.yaml'
publish_dir: s3://openproblems-data/resources/datasets/multimodal
HERE

cat > /tmp/nextflow.config << HERE
process {
  withName:'.*publishStatesProc' {
    memory = '16GB'
    disk = '100GB'
  }
}
HERE

tw launch https://github.com/openproblems-bio/datasets.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/multimodal/load_openproblems_neurips2022_pbmc/main.nf \
  --workspace 53907369739130 \
  --compute-env 1pK56PjjzeraOOC2LDZvN2 \
  --params-file "$params_file" \
  --config /tmp/nextflow.config \
  --labels openproblems_neurips2022_pbmc,dataset_loader \
