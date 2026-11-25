# /// script
# requires-python = ">=3.12"
# dependencies = ["lamindb[aws,bionty]==2.12.2", "anndata==0.13.2", "pyyaml"]
# ///
"""Test resources: mouse brain coronal section.

A coronal section of an FFPE mouse brain, profiled with 10x Genomics Visium CytAssist and
the mouse whole transcriptome probe set (10x Genomics, 2022,
https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard). The
section has 2,310 spots under tissue, with a median of 6,736 genes and 24,862 UMIs per
spot, and comes with an H&E image.

Run from the root of the repository:

    uv run scripts/create_test_resources/mouse_brain_coronal.py
"""

import anndata as ad
import lamindb as ln

import helpers

helpers.check_instance()
ln.track("tsSl1e8PA23J", project=helpers.PROJECT)
helpers.link_to_github("scripts/create_test_resources/mouse_brain_coronal.py")

# remove spots with fewer than 200 detected genes, genes detected in fewer than 50 spots and
# mitochondrial genes, subsample to 600 spots and 500 genes, keeping the spot coordinates
# and the tissue image, and normalise the counts to log CP10k
tenx = "https://cf.10xgenomics.com/samples/spatial-exp/2.0.0/CytAssist_FFPE_Mouse_Brain_Rep1"
output_dir = helpers.run_workflow("spatial/process_tenx_visium", {
  "id": "mouse_brain_coronal",
  "input_expression": f"{tenx}/CytAssist_FFPE_Mouse_Brain_Rep1_filtered_feature_bc_matrix.h5",
  "input_spatial": f"{tenx}/CytAssist_FFPE_Mouse_Brain_Rep1_spatial.tar.gz",
  "dataset_name": "10X Visium - Mouse Brain Coronal",
  "dataset_url": "https://www.10xgenomics.com/datasets/mouse-brain-coronal-section-1-ffpe-2-standard",
  "dataset_summary": "Gene expression library of Mouse Brain (CytAssist FFPE) using the Mouse Whole Transcriptome Probe Set",
  "dataset_description": (
    "FFPE Mouse Brain tissue blocks sectioned as described in Visium CytAssist Spatial Gene Expression for "
    "FFPE - Tissue Preparation Guide Demonstrated Protocol. The H&E stained glass slide with tissue section "
    "was processed via Visium CytAssist instrument to transfer analytes to a Visium CytAssist Spatial Gene "
    "Expression slide. The probe extension and library construction steps follow the standard Visium for "
    "FFPE workflow outside of the instrument. The H&E image was acquired using Olympus VS200 Slide Scanning "
    "Microscope. Sequencing depth was 53,497 reads per spot. Sequencing configuration: 28bp read 1 (16bp "
    "Visium spatial barcode, 12bp UMI), 90bp read 2 (transcript), 10bp i7 sample barcode and 10bp i5 sample "
    "barcode. Key metrics include: 2,310 spots detected under tissue; 6,736 median genes per spot; 24,862 "
    "median UMI counts per spot."
  ),
  "dataset_reference": "10x2022brain",
  "dataset_organism": "mus_musculus",
  "spot_filter_min_genes": 200,
  "gene_filter_min_spots": 50,
  "remove_mitochondrial": True,
  "do_subsample": True,
  "n_obs": 600,
  "n_vars": 500,
  "seed": 123,
  "normalization_methods": ["log_cp10k"],
  "output_dataset": "$id/dataset.anndata.zarr",
  "output_state": "$id/state.yaml",
})

print(ad.read_zarr(output_dir / "mouse_brain_coronal" / "dataset.anndata.zarr"))

# store the processed dataset as an artifact, in a collection
helpers.save_collection(
  output_dir / "mouse_brain_coronal",
  ["dataset.anndata.zarr"],
  key="resources_test/datasets/mouse_brain_coronal",
  description="10x Visium mouse brain coronal section, subsampled to 600 spots and 500 genes",
  # the annotation script has no reader for spatial datasets yet
  annotate=False,
)

ln.finish()
