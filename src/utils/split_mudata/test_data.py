"""A small multimodal MuData store for the tests of split_mudata and combine_mudata."""

import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse


def make_mudata(path, mods=("rna", "prot")):
  rng = np.random.default_rng(0)
  cells = [f"cell{i}" for i in range(10)]
  # a protein named like a gene, which must not be renamed
  var_names = {"rna": ["CD4", "GENE1", "GENE2", "GENE3", "GENE4"], "prot": ["CD4", "CD8", "CD19"], "atac": ["chr1:1-100", "chr1:200-300"]}
  adatas = {}
  for name in mods:
    counts = sparse.csr_matrix(rng.poisson(2, (len(cells), len(var_names[name]))) + 1)
    adatas[name] = ad.AnnData(
      layers={"counts": counts},
      obs=pd.DataFrame({f"{name}_qc": rng.random(len(cells))}, index=cells),
      var=pd.DataFrame({"feature_id": var_names[name], "feature_name": var_names[name]}, index=var_names[name]),
      uns={"source": name},
    )
  mdata = md.MuData(adatas)
  mdata.obs = pd.DataFrame({"cell_type": ["t1", "t2"] * 5, "batch": ["b1"] * 5 + ["b2"] * 5}, index=cells)
  mdata.uns = {"dataset_id": "test", "dataset_name": "Test"}
  mdata.write_zarr(path)
  return mdata
