from os import path
import subprocess
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse

# the source files are not public, so build small ones with the same layout
print(">> Creating input files", flush=True)
rng = np.random.default_rng(0)
n_cells = 20
obs = pd.DataFrame(
  {
    "kaggle_dataset": rng.choice(["train", "test_public", "test_private"], n_cells),
    "day": rng.choice([2, 3, 4], n_cells),
    "donor": rng.choice(["13176", "31800"], n_cells),
    "cell_type": rng.choice(["HSC", "EryP"], n_cells),
  },
  index=[f"cell{i}" for i in range(n_cells)],
)
obs["batch"] = obs["day"].astype(str) + "_" + obs["donor"]

def source(var_names):
  counts = sparse.csr_matrix(rng.poisson(1, (n_cells, len(var_names))).astype(np.float32))
  return ad.AnnData(X=counts.copy(), layers={"counts": counts}, obs=obs.copy(), var=pd.DataFrame(index=var_names))

input_mod1 = "cite_rna_merged.h5ad"
input_mod2 = "cite_prot_merged.h5ad"
source([f"ENSG0000000000{i}_GENE{i}" for i in range(10)]).write_h5ad(input_mod1)
source(["CD86", "CD274", "CD270"]).write_h5ad(input_mod2)

output_file = "output.mudata.zarr"

print(">> Running script", flush=True)
out = subprocess.run(
  [
    meta["executable"],
    "--input_mod1", input_mod1,
    "--input_mod2", input_mod2,
    "--mod1", "GEX",
    "--mod2", "ADT",
    "--output", output_file,
    "--dataset_id", "openproblems_neurips2022/pbmc_cite",
    "--dataset_name", "OpenProblems NeurIPS2022 CITE-Seq",
    "--dataset_url", "https://www.kaggle.com/competitions/open-problems-multimodal/data",
    "--dataset_reference", "lance2024predicting",
    "--dataset_summary", "A short summary.",
    "--dataset_description", "A couple of paragraphs worth of text.",
  ],
  stderr=subprocess.STDOUT
)

if out.stdout:
    print(out.stdout, flush=True)

if out.returncode:
    print(f"script: '{out.args}' exited with an error.", flush=True)
    exit(out.returncode)

print(">> Checking whether the output exists", flush=True)
assert path.exists(output_file), "Output file does not exist"

print(">> Read output mudata", flush=True)
mdata = md.read_zarr(output_file)
print(mdata, flush=True)

print(">> Check that the output fits the expected API", flush=True)
assert set(mdata.mod) == {"rna", "prot"}, f"Expected modalities rna and prot, got {list(mdata.mod)}"
for name, mod in mdata.mod.items():
  assert mod.X is None, f".X is not empty in modality {name}"
  assert "counts" in mod.layers, f"'counts' not found in the layers of modality {name}"
  assert "cell_type" not in mod.obs.columns, f"cell_type is duplicated in the obs of modality {name}"
assert list(mdata["rna"].var["feature_id"]) == [f"ENSG0000000000{i}" for i in range(10)]
assert list(mdata["rna"].var["feature_name"]) == [f"GENE{i}" for i in range(10)]
assert list(mdata["prot"].var["feature_name"]) == ["CD86", "CD274", "CD270"]
for col in ["cell_type", "batch", "donor", "day", "is_train"]:
  assert col in mdata.obs.columns, f"{col} column not found in global obs"
assert mdata.obs["organism_ontology_term_id"].unique().tolist() == ["NCBITaxon:9606"], "organism_ontology_term_id is not NCBITaxon:9606"
assert all(":" not in col for col in mdata.obs.columns), "global obs contains prefixed columns"
assert mdata.uns["dataset_id"] == "openproblems_neurips2022/pbmc_cite"
assert list(mdata.uns["dataset_reference"]) == ["lance2024predicting"]
assert "dataset_organism" not in mdata.uns

print(">> All tests passed successfully", flush=True)
