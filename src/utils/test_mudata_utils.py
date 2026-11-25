import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
import pytest

from mudata_utils import MODALITIES, is_mudata, to_mudata


def _adata(obs_names, n_vars=4, **obs_columns):
  return ad.AnnData(
    layers={"counts": np.ones((len(obs_names), n_vars), dtype=np.int32)},
    obs=pd.DataFrame(obs_columns, index=obs_names),
    var=pd.DataFrame(index=[f"f{i}" for i in range(n_vars)]),
  )


CELLS = ["c1", "c2", "c3"]


def test_modalities_maps_source_labels():
  assert MODALITIES == {"GEX": "rna", "ADT": "prot", "ATAC": "atac"}


def test_to_mudata_moves_shared_columns():
  rna = _adata(CELLS, cell_type=["a", "b", "a"], n_genes=[1, 2, 3])
  prot = _adata(CELLS, cell_type=["a", "b", "a"])
  mdata = to_mudata({"rna": rna, "prot": prot}, uns={"dataset_id": "x"})

  assert list(mdata.obs.columns) == ["cell_type"]
  assert list(mdata["rna"].obs.columns) == ["n_genes"]
  assert list(mdata["prot"].obs.columns) == []
  assert mdata.uns["dataset_id"] == "x"


def test_to_mudata_keeps_columns_that_differ():
  rna = _adata(CELLS, quality=[1, 2, 3])
  prot = _adata(CELLS, quality=[3, 2, 1])
  mdata = to_mudata({"rna": rna, "prot": prot}, uns={})

  assert "quality" not in mdata.obs.columns
  assert list(mdata["rna"].obs["quality"]) == [1, 2, 3]


def test_to_mudata_rejects_different_cells():
  rna = _adata(CELLS)
  prot = _adata(["c1", "c2", "c4"])
  with pytest.raises(ValueError, match="Modalities have different cells"):
    to_mudata({"rna": rna, "prot": prot}, uns={})


def test_to_mudata_has_no_prefixed_columns(tmp_path):
  rna = _adata(CELLS, cell_type=["a", "b", "a"], n_genes=[1, 2, 3])
  prot = _adata(CELLS, cell_type=["a", "b", "a"])
  path = tmp_path / "x.mudata.zarr"
  to_mudata({"rna": rna, "prot": prot}, uns={}).write_zarr(path)

  back = md.read_zarr(path)
  assert list(back.obs.columns) == ["cell_type"]
  assert all(":" not in col for col in back.obs.columns)
  assert all(":" not in col for col in back["rna"].obs.columns)


def test_is_mudata(tmp_path):
  rna = _adata(CELLS)
  to_mudata({"rna": rna}, uns={}).write_zarr(tmp_path / "x.mudata.zarr")
  rna.write_zarr(tmp_path / "x.anndata.zarr")

  assert is_mudata(str(tmp_path / "x.mudata.zarr"))
  assert not is_mudata(str(tmp_path / "x.anndata.zarr"))


def test_to_mudata_drops_dataset_level_uns_of_modalities():
  rna = _adata(CELLS)
  rna.uns = {"dataset_id": "source_id", "organism": "human", "genome": "GRCh38"}
  mdata = to_mudata({"rna": rna}, uns={"dataset_id": "x"})

  assert dict(mdata["rna"].uns) == {"genome": "GRCh38"}
  assert mdata.uns["dataset_id"] == "x"
