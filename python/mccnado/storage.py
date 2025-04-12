import h5py
import os
from pathlib import Path
from typing import List, Optional, Union, Callable

class CoolerBinsLinker:
    def __init__(self, file_path):
        self.file_path = file_path
        self.first_bins_ref = None

    def _recursive_link_bins(self, group):
        for key, item in group.items():
            if isinstance(item, h5py.Group):
                # Check if this group is a 'resolutions' group
                if key == 'resolutions':
                    # Iterate through the resolution groups (e.g., '100', '200', etc.)
                    for resolution_key, resolution_group in item.items():
                        if 'bins' in resolution_group:
                            if self.first_bins_ref is None:
                                # This is the first 'bins' dataset we've encountered
                                self.first_bins_ref = resolution_group['bins']
                            else:
                                # Create a hard link to the first 'bins' dataset
                                del resolution_group['bins']
                                resolution_group['bins'] = self.first_bins_ref
                else:
                    # Recursively process other subgroups
                    self._recursive_link_bins(item)

    def link_bins(self):
        with h5py.File(self.file_path, 'r+') as hdf5_file:
            # Start the recursive linking process from the root group
            self._recursive_link_bins(hdf5_file)


class CoolerMerger:
    def __init__(self, source_paths: List[Union[str, Path]], target_path: Union[str, Path], group_namer: Optional[Callable[[str], str]] = None):
        """
        :param source_paths: List of source HDF5 file paths.
        :param target_path: Path to the target HDF5 file.
        :param group_namer: Optional function(source_path) -> group_name in target file.
                            If None, uses the filename (without extension) as the group name.
        """
        self.source_paths = source_paths
        self.target_path = target_path
        self.group_namer = group_namer or self._default_group_namer

    def _default_group_namer(self, source_path) -> str:
        return Path(source_path).stem

    def merge(self):
        with h5py.File(self.target_path, 'a') as tgt:
            for source_path in self.source_paths:
                group_name = self.group_namer(source_path)
                print(f"Merging '{source_path}' into group '{group_name}'...")
                with h5py.File(source_path, 'r') as src:
                    if group_name in tgt:
                        tgt_group = tgt[group_name]
                    else:
                        tgt_group = tgt.create_group(group_name)
                    self._copy_contents(src, tgt_group)

    def _copy_contents(self, src_group, tgt_group):
        self._copy_attributes(src_group, tgt_group)

        for name, item in src_group.items():
            if isinstance(item, h5py.Dataset):
                if name in tgt_group:
                    del tgt_group[name]  # Overwrite existing
                tgt_group.copy(item, name)
                self._copy_attributes(item, tgt_group[name])
            elif isinstance(item, h5py.Group):
                if name in tgt_group:
                    sub_tgt_group = tgt_group[name]
                else:
                    sub_tgt_group = tgt_group.create_group(name)
                self._copy_contents(item, sub_tgt_group)

    def _copy_attributes(self, src_obj, tgt_obj):
        for key, value in src_obj.attrs.items():
            tgt_obj.attrs[key] = value

