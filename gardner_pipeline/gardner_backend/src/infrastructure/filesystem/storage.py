import json
import pickle
import shutil
from pathlib import Path
from typing import Union, Any, Optional, List

from anndata import AnnData
from lotus.io import read_h5ad
import pandas as pd
from ..workspace_context import workspace_path_manager
from .constants.filesystem_constants import USER_PROJECT_ROOT, RAW_DATA
from loguru import logger


class AssetStorage:
    """
    Handles physical file I/O operations for biological assets.
    The Service Layer uses this class to save/load data without managing absolute paths directly.
    """

    def save_raw_dataset(self, adata: AnnData, project_id: str, dataset_id: str, file_name: str) -> str:
        """
        Saves the initial raw dataset into the 'raw_data' subfolder.
        Structure: PROJECT_ROOT / project_id / dataset_id / raw_data / file_name

        Args:
            adata (AnnData): The object to save.
            project_id (str): Parent project ID.
            dataset_id (str): The specific dataset ID.
            file_name (str): The filename (e.g., 'ds_xxx.h5ad').

        Returns:
            str: The standardized POSIX relative path key for database storage.
        """
        # Construct path: projects/p_xxx/ds_xxx/raw_data/filename.h5ad
        relative_path_obj = Path(USER_PROJECT_ROOT) / project_id / dataset_id / RAW_DATA / file_name

        # Convert to POSIX style (forward slashes)
        relative_key = relative_path_obj.as_posix()

        # Use the generic save method
        return self.save_anndata(adata, relative_key)

    def save_anndata(self, adata: AnnData, relative_key: str) -> str:
        """
        Saves an AnnData object to the file system using a relative key.

        Args:
            adata (AnnData): The object to save.
            relative_key (str): The intended logical path.

        Returns:
            str: The standardized POSIX relative path string.
        """
        # 1. Resolve the absolute path using the Context Manager
        target_abs_path = workspace_path_manager.resolve(relative_key)

        # 2. Ensure parent directory exists
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        # 3. Write the file (using gzip to save disk space)
        logger.info(f"[IO] Saving AnnData to: {target_abs_path}")
        try:
            adata.write_h5ad(target_abs_path, compression='gzip')
        except Exception as e:
            # Clean up if write fails to avoid corrupt files
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to write AnnData file: {e}")

        # 4. Return POSIX path
        return Path(relative_key).as_posix()

    def load_anndata(self, relative_key: str) -> AnnData:
        """
        Loads an AnnData object from the workspace.

        Args:
            relative_key (str): The relative path stored in SQLite.

        Returns:
            AnnData: The loaded data object.
        """
        source_abs_path = workspace_path_manager.resolve(relative_key)

        if not source_abs_path.exists():
            raise FileNotFoundError(f"Asset file missing at: {source_abs_path}")

        return read_h5ad(source_abs_path)

    def save_file(self, content: Union[str, bytes, dict, Any], relative_key: str, **kwargs) -> str:
        """
        Generic method to save arbitrary files to the workspace.
        """
        target_abs_path = workspace_path_manager.resolve(relative_key)
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            # 1. Matplotlib Figure Support
            if hasattr(content, 'savefig'):
                save_kwargs = {'bbox_inches': 'tight'}
                dpi = kwargs.pop('dpi', 600)
                save_kwargs.update(kwargs)
                content.savefig(target_abs_path, dpi=dpi, **save_kwargs)

            # 2. PIL Image Support
            elif hasattr(content, 'save') and not isinstance(content, (bytes, str)):
                content.save(target_abs_path, **kwargs)

            # 3. JSON Dictionary
            elif isinstance(content, (dict, list)) and target_abs_path.suffix == '.json':
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    indent = kwargs.get('indent', 2)
                    json.dump(content, f, indent=indent)

            # 4. Binary Data (Bytes)
            elif isinstance(content, bytes):
                with open(target_abs_path, 'wb') as f:
                    f.write(content)

            # 5. Text Data (String)
            elif isinstance(content, str):
                with open(target_abs_path, 'w', encoding='utf-8') as f:
                    f.write(content)

            # 6. Fallback: Python Pickle
            else:
                with open(target_abs_path, 'wb') as f:
                    pickle.dump(content, f)

        except Exception as e:
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to save generic file ({relative_key}): {e}")

        return Path(relative_key).as_posix()

    def save_incremental_anndata(
            self,
            adata_source: AnnData,
            obs_cols: List[str],
            var_cols: List[str],
            relative_key: str
    ) -> Optional[str]:
        """
        Creates and saves a lightweight AnnData object containing only specific columns.
        """
        try:
            adata_cache = AnnData(
                X=None,
                obs=adata_source.obs[obs_cols].copy(),
                var=adata_source.var[var_cols].copy()
            )

            return self.save_anndata(adata_cache, relative_key)

        except Exception as e:
            logger.warning(f"[IO] Failed to save incremental cache: {e}")
            return None

    def load_and_merge_anndata(self, adata_target: AnnData, relative_key: str) -> bool:
        """
        Loads an incremental cache file and merges its obs/var columns into the target AnnData.
        """
        try:

            try:
                adata_cache = self.load_anndata(relative_key)
            except FileNotFoundError:
                return False

            if not adata_target.obs_names.equals(adata_cache.obs_names):
                logger.warning(f"[IO] Cache index mismatch (Obs) for {relative_key}. Aborting merge.")
                return False

            new_obs_cols = adata_cache.obs.columns.difference(adata_target.obs.columns)
            new_var_cols = adata_cache.var.columns.difference(adata_target.var.columns)

            if len(new_obs_cols) > 0:
                for col in new_obs_cols:
                    adata_target.obs[col] = adata_cache.obs[col]

            if len(new_var_cols) > 0:
                for col in new_var_cols:
                    adata_target.var[col] = adata_cache.var[col]

            logger.info(
                f"[IO] Incremental Merge Success! Added {len(new_obs_cols)} obs cols, {len(new_var_cols)} var cols.")
            return True

        except Exception as e:
            logger.warning(f"[IO] Error during incremental merge: {e}")
            return False

    def save_df_to_csv(self, df: pd.DataFrame, relative_key: str, index: bool = False, **kwargs) -> str:
        """
        Saves a Pandas DataFrame to a CSV file in the workspace.
        """
        target_abs_path = workspace_path_manager.resolve(relative_key)
        target_abs_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            df.to_csv(target_abs_path, index=index, encoding='utf-8', **kwargs)
        except Exception as e:
            if target_abs_path.exists():
                target_abs_path.unlink()
            raise IOError(f"Failed to save CSV file ({relative_key}): {e}")

        return Path(relative_key).as_posix()

    def delete_path(self, relative_key: str) -> bool:
        """
        Deletes a file or a directory (recursively) based on a relative path.

        Args:
            relative_key (str): The relative path to delete.

        Returns:
            bool: True if deletion was successful, False otherwise.
        """
        # Resolve the absolute path
        target_abs_path = workspace_path_manager.resolve(relative_key)

        # Check if path exists
        if not target_abs_path.exists():
            logger.warning(f"[IO] Delete failed: Path does not exist at {target_abs_path}")
            return False

        try:
            if target_abs_path.is_dir():
                logger.info(f"[IO] Deleting directory recursively: {target_abs_path}")
                shutil.rmtree(target_abs_path)
            elif target_abs_path.is_file():
                logger.info(f"[IO] Deleting file: {target_abs_path}")
                target_abs_path.unlink()  # Path.unlink() 用于删除文件
            return True
        except Exception as e:
            logger.error(f"[IO] Failed to delete {target_abs_path}: {e}")
            return False
