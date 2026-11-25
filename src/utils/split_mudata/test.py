import subprocess
import sys
from os import path

import anndata as ad

## VIASH START
meta = {
  "executable": "target/executable/utils/split_mudata/split_mudata",
  "resources_dir": "src/utils/split_mudata",
}
## VIASH END

sys.path.append(meta["resources_dir"])
from test_data import make_mudata

print(">> Split a CITE-seq-like store", flush=True)
raw = make_mudata("raw.mudata.zarr", mods=("rna", "prot"))
subprocess.run(
  [meta["executable"], "--input", "raw.mudata.zarr", "--output_rna", "rna.zarr", "--output_prot", "prot.zarr", "--output_atac", "atac.zarr"],
  check=True,
)

assert path.exists("rna.zarr") and path.exists("prot.zarr"), "rna and prot are written"
assert not path.exists("atac.zarr"), "a modality the store doesn't have is not written"

rna = ad.read_zarr("rna.zarr")
assert list(rna.obs.columns) == ["rna_qc", "cell_type", "batch"], f"global obs is joined in, got {list(rna.obs.columns)}"
assert list(rna.obs["cell_type"]) == list(raw.obs["cell_type"])
assert rna.uns["dataset_id"] == "test" and rna.uns["source"] == "rna", "global uns is copied in"
assert list(rna.var_names) == list(raw["rna"].var_names)
assert "counts" in rna.layers

prot = ad.read_zarr("prot.zarr")
assert list(prot.var_names) == ["CD4", "CD8", "CD19"], "feature names are kept as they are"

print(">> All tests passed successfully", flush=True)
