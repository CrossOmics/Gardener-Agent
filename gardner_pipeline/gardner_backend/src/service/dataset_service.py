from datetime import datetime
from pathlib import Path
from typing import Optional, Union

from fastapi import Depends
from loguru import logger

from infrastructure.filesystem.storage import AssetStorage

from infrastructure.database.model.project_meta_model import ProjectMeta
from infrastructure.database.model.dataset_model import Dataset
from infrastructure.database.dao.dataset_dao import DatasetDAO

from lotus.io import standardize_load
from util.id_generate_utils import generate_filename, generate_business_id


class DatasetService:
    """
    Service layer for managing Dataset lifecycle (Ingestion, Persistence, and DB Records).
    """

    def __init__(self, storage: AssetStorage = Depends(),
                 dataset_dao: DatasetDAO = Depends()):
        self.storage = storage
        self.dataset_dao = dataset_dao

    def import_dataset_from_local(
            self,
            local_file_path: str,
            project: Union[ProjectMeta, str],
            dataset_name: Optional[str] = None,
            organism: Optional[str] = None,  # [New]
            tissue_type: Optional[str] = None,  # [New]
            details: Optional[str] = None  # [New]
    ) -> Optional[Dataset]:
        """
        Orchestrates the user upload process:
        1. Validates & Converts raw data to AnnData (Memory).
        2. Defines a secure relative path in the workspace: {project_id}/{dataset_id}/{filename}.
        3. Saves the standardized .h5ad file.
        4. Inserts the record into the Dataset table via DAO.

        Args:
            local_file_path (str): The absolute path chosen by the user.
            project (ProjectMeta | str): The project object or its ID string.
            dataset_name (str, optional): User-defined name.
            organism (str, optional): Species info.
            tissue_type (str, optional): Tissue info.
            details (str, optional): Detailed description.

        Returns:
            Dataset: The created database model instance, or None if failed.
        """

        # 1. Validate Local File
        local_path_obj = Path(local_file_path)
        if not local_path_obj.exists():
            logger.error(f"[Service] Local file not found: {local_file_path}")
            raise FileNotFoundError(f"File not found: {local_file_path}")

        logger.info(f"Loading AnnData from {local_path_obj}...")
        # 2. Ingest & Convert to Standard AnnData (in memory)
        try:
            adata = standardize_load(str(local_path_obj))

            # Optionally store metadata in AnnData.uns for portability
            if organism:
                adata.uns['organism'] = organism
            if tissue_type:
                adata.uns['tissue_type'] = tissue_type

        except Exception as e:
            logger.error(f"[Service] Failed to load or standardize file: {e}")
            raise ValueError(f"Invalid data format: {e}")

        # 3. Prepare Metadata & Naming Conventions
        original_name = local_path_obj.stem
        dataset_business_id = generate_business_id('ds')

        # Naming convention: {dataset_business_id}_{original_name}.h5ad
        filename = f'{dataset_business_id}_{original_name}.h5ad'
        final_dataset_name = dataset_name if dataset_name else original_name

        # Resolve Project ID string
        if isinstance(project, ProjectMeta):
            project_str_id = project.project_id
        else:
            project_str_id = project

        # 4. Persist File to Workspace (New Hierarchy: Project -> Dataset -> File)
        try:
            # Using the new storage method that utilizes dataset_id as the folder
            saved_relative_path = self.storage.save_raw_dataset(
                adata=adata,
                project_id=project_str_id,
                dataset_id=dataset_business_id,  # This creates the specific folder
                file_name=filename
            )
        except Exception as e:
            logger.error(f"[Service] Storage write failed: {e}")
            raise e

        logger.info(f"[Service] Inserting dataset record into DB: {final_dataset_name}")

        # 5. Insert Record to DB
        dataset_record = DatasetDAO.create_dataset(
            project=project_str_id,
            dataset_id=dataset_business_id,
            dataset_name=final_dataset_name,
            dataset_path=saved_relative_path,
            organism=organism,
            tissue_type=tissue_type,
            details=details
        )

        if dataset_record:
            logger.success(f"[Service] Dataset imported successfully. ID: {dataset_record.id}")
            return dataset_record
        else:
            logger.error("[Service] Database insertion failed.")
            return None

    def get_dataset_by_id(self, pk_id: int) -> Optional[Dataset]:
        """Retrieves an ACTIVE dataset by PK."""
        return DatasetDAO.get_dataset_by_id(pk_id)


    def update_dataset_name(self, dataset_id: str, new_name: str) -> bool:
        """
        Updates the display name of a specific dataset.

        Args:
            dataset_id: The unique business identifier of the dataset.
            new_name: The new display name.

        Returns:
            bool: True if the update was successful.
        """
        # Verify dataset exists
        dataset = self.dataset_dao.get_dataset_by_business_id(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found.")

        # update dataset name
        success = self.dataset_dao.update_dataset(
            pk_id=dataset.id,
            dataset_name=new_name
        )

        if not success:
            logger.error(f"Failed to update dataset name for {dataset_id}")

        return success

    def get_dataset_name(self, dataset_id: str) -> str:
        """
        Retrieves the display name of a specific dataset.
        """
        dataset = self.dataset_dao.get_dataset_by_business_id(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found.")

        return dataset.dataset_name
