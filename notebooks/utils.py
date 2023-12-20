from pathlib import Path
import cellxgene_census
import anndata as ad


def download_data(
    dataset_ids: list,
    dataset_name: str,
    data_dir: str = 'data',
    overwrite: bool = False
) -> ad.AnnData:
    """
    Download data from cellxgene-census and save to h5ad file.
    :param dataset_ids: list of dataset ids
    :param dataset_name: name of the dataset that you want to save the file to
    :param data_dir: directory to save the data
    :param overwrite: whether to overwrite the existing file
    """
    output_file = f'{data_dir}/{dataset_name}.h5ad'
    if overwrite and Path(output_file).exists():
        Path(output_file).unlink()
    elif Path(output_file).exists():
        return ad.read_h5ad(output_file)
    
    if len(dataset_ids) > 1:
        adatas = []
        for dataset_id in dataset_ids:
            per_id_file = f'{data_dir}/{dataset_name}_{dataset_id}.h5ad'
            print(f'download to {per_id_file}...')
            cellxgene_census.download_source_h5ad(dataset_id, to_path=per_id_file)
            adatas.append(ad.read_h5ad(per_id_file))
        print('concatenate...')
        adata = ad.concat(adatas)
        print(f'write to {output_file}...')
        adata.write(output_file)
        return adata
    else:
        print(f'download to {output_file}...')
        cellxgene_census.download_source_h5ad(dataset_ids[0], to_path=output_file)
    adata = ad.read_h5ad(output_file)
    return adata