"""Tests for annotate_artifacts.py that don't need a LaminDB instance.

    uv run --with pytest --with "lamindb[aws,bionty]==2.12.2" --with anndata==0.13.2 \
      --with mudata==0.4.1 --with pyyaml pytest scripts/lamin/test_annotate_artifacts.py
"""

import sys
from pathlib import Path

import bionty as bt

sys.path.insert(0, str(Path(__file__).parent))
import annotate_artifacts  # noqa: E402


def test_placeholder_registries_cover_cellxgene_unknowns():
  # CELLxGENE schema 5 allows "unknown" in these columns (and "na" for ethnicity), so they
  # need a placeholder term to validate; the Census release of 2025-11-08 has unknown cell types
  assert annotate_artifacts.PLACEHOLDER_REGISTRIES == {
    "cell_type_ontology_term_id": bt.CellType,
    "development_stage_ontology_term_id": bt.DevelopmentalStage,
    "self_reported_ethnicity_ontology_term_id": bt.Ethnicity,
    "sex_ontology_term_id": bt.Phenotype,
  }
