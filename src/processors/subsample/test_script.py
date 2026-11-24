import sys
import os
import pytest
import anndata as ad
import mudata as md
import numpy as np
import pandas as pd
from scipy import sparse

## VIASH START
meta = {
    "resources_dir": "resources_test/common"
}
## VIASH END

input_path = f"{meta['resources_dir']}/cxg_mouse_pancreas_atlas/dataset.anndata.zarr"
input = ad.read_zarr(input_path)

def test_even_sampling(run_component):
    output_path = "output.zarr"
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--even",
        "--seed", "123",
        "--n_obs", "100",
        "--n_vars", "120"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_zarr(output_path)

    assert output.n_obs <= 100, "n_obs should be <= 100"
    assert output.n_vars <= 120, "n_vars should be <= 100"


def test_keep_functionality(run_component):
    output_path = "output.zarr"

    # take categories and well-expressed genes from the test data
    cell_types = input.obs["cell_type"].value_counts().index[:2].tolist()
    batches = input.obs["batch"].value_counts().index[:2].tolist()
    expression = np.asarray(input.layers["counts"].sum(axis=0)).ravel()
    keep_features = input.var_names[np.argsort(-expression)[:3]].tolist()

    run_component([
        "--input", input_path,
        "--keep_cell_type_categories", ";".join(cell_types),
        "--keep_batch_categories", ";".join(batches),
        "--keep_features", ";".join(keep_features),
        "--output", output_path,
        "--seed", "123"
    ])

    # Checking whether file exists
    assert os.path.exists(output_path), "Output file not found"

    # Check that test output fits expected API
    output = ad.read_zarr(output_path)

    assert output.n_obs <= 500, "n_obs should be <= 500"
    assert output.n_vars <= 500, "n_vars should be <= 500"
    for feat in keep_features:
        assert feat in output.var_names, f"{feat} should be in output.var_names"
    assert set(output.obs["cell_type"]) <= set(cell_types), "only the kept cell types should remain"
    assert set(output.obs["batch"]) <= set(batches), "only the kept batches should remain"


def _mudata_input(path, n_cells=200, prot_cells=None):
    rng = np.random.default_rng(0)
    cells = [f"cell{i}" for i in range(n_cells)]

    def modality(prefix, n_vars, obs_names):
        counts = sparse.csr_matrix(rng.poisson(1, (len(obs_names), n_vars)) + 1)
        return ad.AnnData(
            layers={"counts": counts},
            obs=pd.DataFrame(index=obs_names),
            var=pd.DataFrame(index=[f"{prefix}{i}" for i in range(n_vars)]),
        )

    mdata = md.MuData({
        "rna": modality("gene", 3000, cells),
        "prot": modality("protein", 140, prot_cells or cells),
    })
    mdata.obs["batch"] = pd.Categorical(np.resize(["b1", "b2", "b3", "b4"], mdata.n_obs))
    mdata.obs["cell_type"] = "t1"
    mdata.uns["dataset_id"] = "test"
    mdata.write_zarr(path)
    return path


def test_subsample_mudata_per_modality(run_component):
    input_path = _mudata_input("input.mudata.zarr")
    outputs = []
    for output_path in ["output1.mudata.zarr", "output2.mudata.zarr"]:
        run_component([
            "--input", input_path,
            "--output", output_path,
            "--even",
            "--seed", "1",
            "--n_obs", "100",
            "--n_vars", "1500",
        ])
        outputs.append(md.read_zarr(output_path))

    output = outputs[0]
    assert set(output.mod) == {"rna", "prot"}
    assert output.n_obs == 100
    assert output["rna"].n_vars == 1500, "n_vars applies to each modality"
    assert output["prot"].n_vars == 140, "a small modality keeps all its features"
    assert list(output["rna"].obs_names) == list(output["prot"].obs_names)
    assert output.uns["dataset_id"] == "test_subsample"
    assert "batch" in output.obs and "batch" not in output["rna"].obs
    assert list(output.obs_names) == list(outputs[1].obs_names), "same seed, same cells"
    assert list(output["rna"].var_names) == list(outputs[1]["rna"].var_names), "same seed, same features"


def test_subsample_mudata_keep_features(run_component):
    input_path = _mudata_input("input.mudata.zarr")
    output_path = "output.mudata.zarr"
    run_component([
        "--input", input_path,
        "--output", output_path,
        "--seed", "1",
        "--n_obs", "50",
        "--n_vars", "20",
        "--keep_features", "gene2999;protein7",
    ])
    output = md.read_zarr(output_path)

    assert "gene2999" in output["rna"].var_names
    assert "protein7" in output["prot"].var_names
    assert output["rna"].n_vars == 20 and output["prot"].n_vars == 20


def test_subsample_mudata_rejects_different_cells(run_component):
    cells = [f"cell{i}" for i in range(200)]
    input_path = _mudata_input("input.mudata.zarr", prot_cells=cells[:100] + [f"other{i}" for i in range(100)])
    with pytest.raises(Exception):
        run_component([
            "--input", input_path,
            "--output", "output.mudata.zarr",
            "--seed", "1",
        ])

if __name__ == '__main__':
    sys.exit(pytest.main([__file__, "--capture=no"], plugins=["viashpy"]))
