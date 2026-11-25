"""Tests for annotate_artifacts.py that don't need a LaminDB instance.

    uv run --with pytest --with "lamindb[aws,bionty]==2.12.2" --with anndata==0.13.2 \
      --with mudata==0.4.1 --with pyyaml pytest scripts/lamin/test_annotate_artifacts.py
"""

import sys
from pathlib import Path

import bionty as bt

sys.path.insert(0, str(Path(__file__).parent))
import annotate_artifacts  # noqa: E402


def test_describe_dataframe_with_nullable_string_index(tmp_path):
  # with pandas 3, anndata writes string indexes as nullable string arrays (a group), as
  # the OP3 loader does
  import anndata as ad
  import numpy as np
  import pandas as pd
  import zarr

  obs = pd.DataFrame({"a": [1, 2, 3], "b": ["x", "y", "z"]}, index=pd.Index(["c1", "c2", None], dtype="string"))
  with ad.settings.override(allow_write_nullable_strings=True):
    ad.AnnData(X=np.zeros((3, 2)), obs=obs).write_zarr(tmp_path / "test.zarr")
  root = zarr.open_group(tmp_path / "test.zarr", mode="r")
  assert isinstance(root["obs"][root["obs"].attrs["_index"]], zarr.Group)

  assert annotate_artifacts._describe_element(root["obs"]) == {"type": "dataframe", "shape": [3, 2]}


def test_placeholder_registries_cover_cellxgene_unknowns():
  # CELLxGENE schema 5 allows "unknown" in these columns (and "na" for ethnicity), so they
  # need a placeholder term to validate; the Census release of 2025-11-08 has unknown cell types
  assert annotate_artifacts.PLACEHOLDER_REGISTRIES == {
    "cell_type_ontology_term_id": bt.CellType,
    "development_stage_ontology_term_id": bt.DevelopmentalStage,
    "self_reported_ethnicity_ontology_term_id": bt.Ethnicity,
    "sex_ontology_term_id": bt.Phenotype,
  }
