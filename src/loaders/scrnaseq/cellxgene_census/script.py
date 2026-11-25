import cellxgene_census
import pandas as pd
import scanpy as sc
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

## VIASH START
par = {
    "input_uri": None,
    "census_version": "stable",
    "species": "mus_musculus",
    "obs_value_filter": "dataset_id == '49e4ffcc-5444-406d-bdee-577127404ba8'",
    "cell_filter_grouping": None,
    "cell_filter_minimum_count": None,
    "obs_batch": [ "donor_id" ],
    "obs_batch_separator": "+",
    "dataset_name": "pretty name",
    "dataset_url": "url",
    "dataset_reference": "ref",
    "dataset_summary": "summ",
    "dataset_description": "desc",
    "output": "output.zarr",
}
## VIASH END

# organisms in the Census, with their NCBI Taxonomy term
CENSUS_ORGANISMS = {
    "homo_sapiens": ("Homo sapiens", "NCBITaxon:9606"),
    "mus_musculus": ("Mus musculus", "NCBITaxon:10090"),
}

def connect_census(uri, census_version):
    """
    Connect to CellxGene Census or user-provided TileDBSoma object
    """
    ver = census_version or "stable"
    logger.info("Connecting to CellxGene Census at %s", f"'{uri}'" if uri else f"version '{ver}'")
    return cellxgene_census.open_soma(uri=uri, census_version=ver)

def get_anndata(census_connection, par):
    logger.info("Getting gene expression data based on `%s` query.", par["obs_value_filter"])
    return cellxgene_census.get_anndata(
        census=census_connection,
        obs_value_filter=par["obs_value_filter"],
        organism=par["species"]
    )

def filter_min_cells_per_group(adata, par):
    n_cells_before, _ = adata.shape
    cell_count = adata.obs \
        .groupby(par["cell_filter_grouping"])["soma_joinid"] \
        .transform("count")
    adata = adata[cell_count >= par["cell_filter_minimum_count"]]
    n_cells_after, _ = adata.shape
    logger.info(
        "Removed %s cells based on %s cell_filter_minimum_count of %s cell_filter_grouping."
        % ((n_cells_before - n_cells_after), par["cell_filter_minimum_count"], par["cell_filter_grouping"])
    )
    return adata

def filter_by_counts(adata, par):
    logger.info("Remove cells with few counts and genes with few counts.")
    n_cells_before, n_genes_before = adata.shape
    # remove cells with few counts and genes with few counts
    scanpy_proc = {
        par["cell_filter_min_counts"]: (sc.pp.filter_cells, "min_counts"),
        par["cell_filter_min_genes"]: (sc.pp.filter_cells, "min_genes"),
        par["gene_filter_min_counts"]: (sc.pp.filter_genes, "min_counts"),
        par["gene_filter_min_cells"]: (sc.pp.filter_genes, "min_cells"),
    }
    for threshold, (func, arg) in scanpy_proc.items():
        if threshold:
            func(adata, **{arg: threshold})
    n_cells_after, n_genes_after = adata.shape
    logger.info("Removed %s cells and %s genes.", (n_cells_before - n_cells_after), (n_genes_before - n_genes_after))

def move_x_to_layers(adata):
    logger.info("Move .X to .layers['counts']")
    adata.layers["counts"] = adata.X
    adata.X = None

def add_batch_to_obs(adata, par):
    logger.info("Add batch to the AnnData object.")
    if par["obs_batch"]:
        # fetch batch columns from obs
        cols = [adata.obs[key] for key in par["obs_batch"]]
        
        # join cols
        obs_batch = [par["obs_batch_separator"].join(row) for row in zip(*cols)]

        # store in adata
        adata.obs["batch"] = obs_batch

def add_organism_to_obs(adata, par):
    logger.info("Add organism to obs.")
    species = par["species"].lower().replace(" ", "_")
    if species not in CENSUS_ORGANISMS:
        raise ValueError(f"Unknown species '{par['species']}', expected one of {list(CENSUS_ORGANISMS)}")
    label, term_id = CENSUS_ORGANISMS[species]
    adata.obs["organism"] = pd.Categorical([label] * adata.n_obs)
    adata.obs["organism_ontology_term_id"] = pd.Categorical([term_id] * adata.n_obs)

def add_metadata_to_uns(adata, par):
    logger.info("Add metadata to the AnnData object.")
    for key in ["dataset_id", "dataset_name", "dataset_url", "dataset_reference", "dataset_summary", "dataset_description"]:
        adata.uns[key] = par[key]

def print_unique(adata, column):
    formatted = "', '".join(adata.obs[column].unique())
    logger.info(f"Unique {column}: ['{formatted}']")

def print_summary(adata):
    logger.info(f"Resulting dataset: {adata}")

    logger.info("Summary of dataset:")
    obs_fields = ["assay", "assay_ontology_term_id", "cell_type", "cell_type_ontology_term_id", "dataset_id", "development_stage", "development_stage_ontology_term_id", "disease", "disease_ontology_term_id", "organism", "organism_ontology_term_id", "tissue", "tissue_ontology_term_id", "tissue_general", "tissue_general_ontology_term_id"]
    for field in obs_fields:
        print_unique(adata, field)
def write_anndata(adata, par):
    logger.info("Writing AnnData object to '%s'", par["output"])

    adata.write_zarr(par["output"])

def main(par, meta):
    # check arguments
    if (par["cell_filter_grouping"] is None) != (par["cell_filter_minimum_count"] is None):
        raise NotImplementedError(
            "You need to specify either both or none of the following parameters: cell_filter_grouping, cell_filter_minimum_count"
        )
    
    with connect_census(uri=par["input_uri"], census_version=par["census_version"]) as conn:
        adata = get_anndata(conn, par)
    
    print(f"AnnData: {adata}", flush=True)

    if par["cell_filter_grouping"] is not None:
        adata = filter_min_cells_per_group(adata, par)

    # remove cells with few counts and genes with few counts
    filter_by_counts(adata, par)

    # logger.log(f"Filtered AnnData: {adata}")
    print(f"Filtered AnnData: {adata}", flush=True)

    # use feature_id as var_names
    adata.var_names = adata.var["feature_id"]
    adata.var.index.name = None

    # not needed as long as we have our own implementation of `get_anndata`
    # move .X to .layers["counts"]
    move_x_to_layers(adata)

    # add batch to obs
    add_batch_to_obs(adata, par)

    # add organism to obs
    add_organism_to_obs(adata, par)

    # add metadata to uns
    add_metadata_to_uns(adata, par)

    # print summary
    print_summary(adata)

    # write output to file
    write_anndata(adata, par)


if __name__ == "__main__":
    main(par, meta)
