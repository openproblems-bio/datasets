"""Tests for process_datasets.py, without a LaminDB instance.

    uv run --with pytest --with pyyaml pytest scripts/create_resources
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import process_datasets  # noqa: E402

UID = "abcdefghijklmnop0000"


def test_scrnaseq_key():
  assert process_datasets.raw_key_to_param("scrnaseq/cellxgene_census/hlca/raw.anndata.zarr", UID, "openproblems/data-dev") == {
    "id": "cellxgene_census/hlca",
    "input": f"lamin://openproblems/data-dev/artifact/{UID}",
    "output_dataset": "$id/dataset.anndata.zarr",
    "output_state": "$id/process_state.yaml",
  }


def test_mudata_key():
  param = process_datasets.raw_key_to_param("multimodal/openproblems_neurips2021/bmmc_cite/raw.mudata.zarr", UID, "openproblems/data-dev")
  assert param["id"] == "openproblems_neurips2021/bmmc_cite"
  assert param["output_dataset"] == "$id/dataset.mudata.zarr"


def test_op3_key():
  param = process_datasets.raw_key_to_param("scrnaseq/openproblems_op3/raw.anndata.zarr", UID, "openproblems/data-dev")
  assert param["id"] == "openproblems_op3"


def test_rejects_processed_key():
  with pytest.raises(ValueError):
    process_datasets.raw_key_to_param("scrnaseq/cellxgene_census/hlca/dataset.anndata.zarr", UID, "openproblems/data-dev")


def test_rejects_test_resource_key():
  with pytest.raises(ValueError):
    process_datasets.raw_key_to_param("resources_test/datasets/cxg_mouse_pancreas_atlas/raw.anndata.zarr", UID, "openproblems/data-dev")


def test_large_ids_come_from_the_load_scripts():
  assert process_datasets.LARGE_IDS == {
    "cellxgene_census/tabula_sapiens",
    "cellxgene_census/human_retina_cell_atlas",
    "openproblems_op3",
    "openproblems_neurips2022/pbmc_multiome",
  }
