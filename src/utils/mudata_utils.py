"""Helpers for the components that read or write multimodal MuData stores."""

from pathlib import Path

import mudata as md

# modality labels used by the sources, and the name of that modality in a MuData store
MODALITIES = {"GEX": "rna", "ADT": "prot", "ATAC": "atac"}


def to_mudata(mods, uns):
  """Combine AnnData objects that measure the same cells into one MuData object.

  Columns of `obs` that are identical in every modality describe the cell rather than
  the measurement, so they are moved to the global `obs`. `uns` becomes the global `uns`.
  """
  names = list(mods)
  first = mods[names[0]].obs
  for name in names[1:]:
    if not mods[name].obs_names.equals(first.index):
      raise ValueError(f"Modalities have different cells: '{names[0]}' and '{name}'")

  shared = [
    col for col in first.columns
    if all(col in mods[name].obs and mods[name].obs[col].equals(first[col]) for name in names[1:])
  ]
  obs = first[shared].copy()
  for adata in mods.values():
    adata.obs = adata.obs.drop(columns=shared)

  mdata = md.MuData(mods)
  mdata.obs = obs
  mdata.uns = dict(uns)
  return mdata


def is_mudata(path):
  """Whether `path` is a MuData zarr store rather than an AnnData one."""
  return (Path(path) / "mod").is_dir()
