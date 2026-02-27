from datetime import datetime
from typing import List, Optional, Dict, Any, Union
from peewee import DoesNotExist, IntegrityError

from ..model.dataset_model import Dataset
from ..model.project_meta_model import ProjectMeta


class DatasetDAO:
    """
    Data Access Object (DAO) for managing Dataset database operations.
    Handles interaction with the 'dataset' table, linking raw files to projects.
    """

    @staticmethod
    def create_dataset(
            project: Union[ProjectMeta, str],
            dataset_id: str,
            dataset_name: str,
            dataset_path: str,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None,
            details: Optional[str] = None,
            ext_info: Optional[Dict[str, Any]] = None
    ) -> Optional[Dataset]:
        """
        Creates a new Dataset record.

        Args:
            project (ProjectMeta | str): The parent project instance or project_id string.
            dataset_id (str): Unique business identifier.
            dataset_name (str): Display name for the dataset.
            dataset_path (str): Relative file path in storage.
            organism (str, optional): Species (e.g., Human).
            tissue_type (str, optional): Tissue (e.g., PBMC).
            details (str, optional): Detailed description.
            ext_info (dict, optional): Dictionary containing file metadata.
        """
        try:
            dataset = Dataset.create(
                project_id=project,
                dataset_id=dataset_id,
                dataset_name=dataset_name,
                dataset_path=dataset_path,
                organism=organism,
                tissue_type=tissue_type,
                details=details,
                uploaded_time=datetime.utcnow(),
                is_deleted=False,
                ext_info=ext_info or {}
            )
            return dataset
        except IntegrityError as e:
            print(f"[Error] Failed to create dataset '{dataset_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating dataset: {e}")
            return None

    @staticmethod
    def get_dataset_by_id(pk_id: int) -> Optional[Dataset]:
        """Retrieves a single ACTIVE dataset by primary key."""
        try:
            return Dataset.get((Dataset.id == pk_id) & (Dataset.is_deleted == False))
        except DoesNotExist:
            print(f"[Warning] Dataset with PK {pk_id} not found.")
            return None

    @staticmethod
    def get_dataset_by_business_id(dataset_id: str) -> Optional[Dataset]:
        """Retrieves an ACTIVE dataset by business ID."""
        try:
            return Dataset.get((Dataset.dataset_id == dataset_id) & (Dataset.is_deleted == False))
        except DoesNotExist:
            print(f"[Warning] Dataset with business ID '{dataset_id}' not found.")
            return None

    @staticmethod
    def get_datasets_by_project(project_pk: str, include_deleted: bool = False) -> List[Dataset]:
        """
        Retrieves datasets linked to a specific project.
        """
        try:
            query = Dataset.select().where(Dataset.project_id == project_pk)
            if not include_deleted:
                query = query.where(Dataset.is_deleted == False)

            return list(query.order_by(Dataset.uploaded_time.desc()))
        except Exception as e:
            print(f"[Error] Failed to retrieve datasets for project {project_pk}: {e}")
            return []

    @staticmethod
    def filter_datasets(
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None
    ) -> List[Dataset]:
        """
        Filters ACTIVE datasets based on organism and/or tissue type.
        """
        try:
            query = Dataset.select().where(Dataset.is_deleted == False)

            if organism:
                query = query.where(Dataset.organism == organism)
            if tissue_type:
                query = query.where(Dataset.tissue_type == tissue_type)

            return list(query.order_by(Dataset.uploaded_time.desc()))
        except Exception as e:
            print(f"[Error] Failed to filter datasets: {e}")
            return []

    @staticmethod
    def update_dataset(
            pk_id: int,
            dataset_name: Optional[str] = None,
            dataset_path: Optional[str] = None,
            organism: Optional[str] = None,
            tissue_type: Optional[str] = None,
            details: Optional[str] = None,
            ext_info_update: Optional[Dict[str, Any]] = None
    ) -> bool:
        """
        Updates an existing dataset record.
        """
        try:
            dataset = Dataset.get_by_id(pk_id)

            if dataset_name is not None:
                dataset.dataset_name = dataset_name
            if dataset_path is not None:
                dataset.dataset_path = dataset_path
            if organism is not None:
                dataset.organism = organism
            if tissue_type is not None:
                dataset.tissue_type = tissue_type
            if details is not None:
                dataset.details = details

            if ext_info_update is not None:
                current_info = dataset.ext_info or {}
                current_info.update(ext_info_update)
                dataset.ext_info = current_info

            dataset.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot update: Dataset ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update dataset {pk_id}: {e}")
            return False

    @staticmethod
    def delete_dataset(pk_id: int, hard_delete: bool = False) -> bool:
        """
        Deletes a dataset record. Defaults to SOFT DELETE.

        Args:
            pk_id (int): Primary key of the dataset.
            hard_delete (bool): If True, physically removes the record.
        """
        try:
            dataset = Dataset.get_by_id(pk_id)
            if hard_delete:
                dataset.delete_instance()
                print(f"[Info] Dataset ID {pk_id} HARD deleted successfully.")
            else:
                dataset.is_deleted = True
                dataset.save()
                print(f"[Info] Dataset ID {pk_id} SOFT deleted successfully.")
            return True
        except DoesNotExist:
            print(f"[Warning] Cannot delete: Dataset ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to delete dataset {pk_id}: {e}")
            return False