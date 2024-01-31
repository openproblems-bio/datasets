import warnings
from pathlib import Path
import cellxgene_census
import anndata as ad
import zarr
import h5py
from scipy.sparse import csr_matrix
from anndata.experimental import read_elem, sparse_dataset


def load_data(
    dataset_ids: list,
    dataset_name: str,
    data_dir: str = 'data',
    overwrite: bool = False,
    **kwargs,
) -> ad.AnnData:
    """
    Download data from cellxgene-census and save to h5ad file, read file and return AnnData.
    If the file already exists and overwrite is False, read the file and return AnnData.
    
    :param dataset_ids: list of dataset ids
    :param dataset_name: name of the dataset that you want to save the file to
    :param data_dir: directory to save the data
    :param overwrite: whether to overwrite the existing file
    :kwargs: additional parameters to pass to read_anndata
    """
    output_file = f'{data_dir}/{dataset_name}.h5ad'
    if overwrite and Path(output_file).exists():
        Path(output_file).unlink()
    elif Path(output_file).exists():
        return read_anndata(output_file, **kwargs)
    
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
    print(f'read {output_file}...')
    return read_anndata(output_file, **kwargs)


def read_anndata(
    file: str,
    backed: bool = False,
    **kwargs
) -> ad.AnnData:
    """
    Read anndata file
    :param file: path to anndata file in zarr or h5ad format
    :param kwargs: AnnData parameter to zarr group mapping
    """
    
    def get_file_reader(file):
        if file.endswith(('.zarr', '.zarr/')):
            func = zarr.open
            file_type = 'zarr'
        elif file.endswith('.h5ad'):
            func = h5py.File
            file_type = 'h5py'
        else:
            raise ValueError(f'Unknown file format: {file}')
        return func, file_type
    
    assert Path(file).exists(), f'File not found: {file}'
    
    func, file_type = get_file_reader(file)
    
    f = func(file, 'r')
    kwargs = {x: x for x in f} if not kwargs else kwargs
    if len(f.keys()) == 0:
        return ad.AnnData()
    # check if keys are available
    for name, slot in kwargs.items():
        if slot not in f:
            warnings.warn(
                f'Cannot find "{slot}" for AnnData parameter `{name}` from "{file}"'
            )
    adata = read_partial(f, backed=backed, **kwargs)
    if not backed and file_type == 'h5py':
        f.close()
    
    return adata


def read_partial(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    force_sparse_types: [str, list] = None,
    **kwargs
) -> ad.AnnData:
    """
    Partially read zarr or h5py groups
    :params group: file group
    :params force_sparse_types: encoding types to convert to sparse_dataset via csr_matrix
    :params backed: read sparse matrix as sparse_dataset
    :params **kwargs: dict of slot_name: slot, by default use all available slot for the zarr file
    :return: AnnData object
    """
    if force_sparse_types is None:
        force_sparse_types = []
    elif isinstance(force_sparse_types, str):
        force_sparse_types = [force_sparse_types]
    slots = {}
    if backed:
        print('Read as backed sparse matrix...')
    
    for slot_name, slot in kwargs.items():
        print(f'Read slot "{slot}", store as "{slot_name}"...')
        if slot not in group:
            warnings.warn(f'Slot "{slot}" not found, skip...')
            slots[slot_name] = None
        else:
            elem = group[slot]
            iospec = ad._io.specs.get_spec(elem)
            if iospec.encoding_type in ("csr_matrix", "csc_matrix") and backed:
                slots[slot_name] = sparse_dataset(elem)
            elif iospec.encoding_type in force_sparse_types:
                slots[slot_name] = csr_matrix(read_elem(elem))
                if backed:
                    slots[slot_name] = sparse_dataset(slots[slot_name])
            else:
                slots[slot_name] = read_elem(elem)
    return ad.AnnData(**slots)

