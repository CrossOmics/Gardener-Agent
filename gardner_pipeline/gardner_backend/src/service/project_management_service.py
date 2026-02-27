from pathlib import Path
from typing import Optional, List

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO
from infrastructure.database.dao.project_meta_dao import ProjectMetaDAO
from infrastructure.database.model.project_meta_model import ProjectMeta
from util.id_generate_utils import generate_business_id
from loguru import logger


class ProjectManagementService:
    """
    Service layer for Project Metadata management.
    Handles business logic such as ID generation, path resolution, and bridging
    data access to the controller layer.
    """

    def __init__(self):
        pass

    def create_project(
            self,
            project_name: str,
            description: Optional[str] = None
    ) -> Optional[ProjectMeta]:
        """
        Creates a new project. Generates a unique Project ID and assigns a workspace path.

        Args:
            project_name: User-defined name.
            description: Optional description.

        Returns:
            ProjectMeta: The created project object, or None if failed.
        """
        # Generate Unique Business ID
        # Format: p_<timestamp>_<short_uuid>
        project_id = generate_business_id('p')

        # Determine Project Root Path
        # Logic: USER_PROJECT_ROOT / project_id
        # We store the relative path string (compatible with pathlib logic later)
        from infrastructure.filesystem.constants.filesystem_constants import USER_PROJECT_ROOT
        project_path_obj = Path(USER_PROJECT_ROOT) / project_id
        project_path_str = project_path_obj.as_posix()

        # Persist to Database via DAO
        # Organism and tissue_type are no longer project-level attributes
        project = ProjectMetaDAO.create_project(
            project_id=project_id,
            project_name=project_name,
            project_path=project_path_str,
            description=description
        )

        if project:
            logger.info(f"Created project: {project_name} ({project_id})")
        else:
            logger.error(f"Failed to create project: {project_name}")

        return project

    def get_project_by_id(self, project_id: str) -> Optional[ProjectMeta]:
        """
        Retrieves an ACTIVE project by its business ID (e.g., 'p_2025...').
        """
        return ProjectMetaDAO.get_project_by_project_id(project_id)

    def get_project_by_pk(self, pk_id: int) -> Optional[ProjectMeta]:
        """
        Retrieves an ACTIVE project by its database primary key.
        """
        return ProjectMetaDAO.get_project_by_id(pk_id)

    def list_all_projects(self, include_deleted: bool = False) -> List[ProjectMeta]:
        """
        Returns a list of projects, ordered by last update.
        """
        return ProjectMetaDAO.get_all_projects(include_deleted=include_deleted)

    def update_project_details(
            self,
            pk_id: int,
            name: Optional[str] = None,
            description: Optional[str] = None,
            ext_info: Optional[str] = None
    ) -> bool:
        """
        Updates project metadata.
        """
        return ProjectMetaDAO.update_project(
            pk_id=pk_id,
            project_name=name,
            description=description,
            ext_info=ext_info
        )

    def delete_project(self, pk_id: int, hard_delete: bool = False) -> bool:
        """
        Deletes a project record.

        Args:
            pk_id: The primary key of the project.
            hard_delete: If True, physically removes the record from DB.
                         If False (default), marks as is_deleted=1.
        """
        return ProjectMetaDAO.delete_project(pk_id, hard_delete=hard_delete)

    def delete_project_by_id(self, project_id: str) -> dict:
        """
        Deletes a project by its business ID.

        Args:
            project_id: Unique string ID (e.g., 'p_123').

        Returns:
            dict: Success message.
        """
        # 1. Retrieve project to ensure it exists and get path
        project = self.get_project_by_id(project_id)
        if not project:
            raise ValueError(f"Project {project_id} not found.")

        # 2. Delete from Database
        db_success = ProjectMetaDAO.delete_project(project.id, hard_delete=True)

        if not db_success:
            raise RuntimeError("Failed to delete project record from database.")

        # 3. Delete from File System
        from infrastructure.filesystem.storage import AssetStorage
        storage = AssetStorage()

        io_success = storage.delete_path(project.project_path)

        if not io_success:
            logger.warning(f"Project record deleted, but failed to cleanup files at: {project.project_path}")

        return {"msg": f"Project {project_id} deleted successfully"}

    def update_project_name(self, project_id: str, new_name: str) -> bool:
        """
        Updates the display name of a specific project based on its business ID.

        Args:
            project_id: The unique business identifier of the project.
            new_name: The new display name.

        Returns:
            bool: True if the update was successful.
        """
        # Verify project exists and get its internal primary key
        project = self.get_project_by_id(project_id)
        if not project:
            raise ValueError(f"Project {project_id} not found.")

        # Call existing DAO update method using the internal PK
        success = self.update_project_details(
            pk_id=project.id,
            name=new_name
        )

        if not success:
            logger.error(f"Failed to update project name for {project_id}")

        return success

    def get_project_name(self, project_id: str) -> str:
        """
        Retrieves the display name of a specific project based on its business ID.
        """
        project = self.get_project_by_id(project_id)
        if not project:
            raise ValueError(f"Project {project_id} not found.")

        return project.project_name

    def get_project_entity_map(self, project_id: str) -> dict:
        """
        Retrieves the full hierarchical map of IDs and Display Names for a project.
        Constructs the requested JSON structure: Project -> Datasets -> Stages -> Snapshots.
        """
        # 1. Get Project
        project = self.get_project_by_id(project_id)
        if not project:
            raise ValueError(f"Project {project_id} not found.")

        result_map = {
            "project_id": project.project_id,
            "project_name": project.project_name,
            "datasets": {}
        }

        datasets = DatasetDAO.get_datasets_by_project(project.project_id)

        for dataset in datasets:
            ds_id = dataset.dataset_id

            # Initialize dataset entry
            result_map["datasets"][ds_id] = {
                "dataset_name": dataset.dataset_name,
                "snapshots": {}
            }

            snapshots = AnalysisSnapshotsDAO.get_snapshots_by_dataset(ds_id)
            for snapshot in snapshots:
                stage = snapshot.branch_name
                snap_id = snapshot.snapshot_id
                snap_name = snapshot.snapshot_name

                # Initialize stage dictionary if it doesn't exist
                if stage not in result_map["datasets"][ds_id]["snapshots"]:
                    result_map["datasets"][ds_id]["snapshots"][stage] = {}

                # Assign snapshot ID -> Name
                result_map["datasets"][ds_id]["snapshots"][stage][snap_id] = snap_name

        return result_map

    def delete_dataset(self, dataset_id: str, keep_final: bool = False) -> dict:
        """
        Deletes a dataset by its business ID.
        If keep_final is True, it retains ONLY the most recent 'Annotation' snapshot
        and deletes all other snapshots. If no 'Annotation' snapshot exists, it
        deletes the entire dataset.

        Args:
            dataset_id: Unique string ID (e.g., 'ds_123').
            keep_final: Flag to retain the latest Annotation snapshot.

        Returns:
            dict: Success message.
        """
        # 1. Retrieve dataset to ensure it exists
        # Note: Ensure DatasetDAO is imported at the top of your service file
        dataset = DatasetDAO.get_dataset_by_business_id(dataset_id)
        if not dataset:
            raise ValueError(f"Dataset {dataset_id} not found.")

        # Determine Dataset File Path
        project_id = dataset.project_id.project_id
        dataset_dir_key = f"projects/{project_id}/{dataset_id}"

        from infrastructure.filesystem.storage import AssetStorage
        storage = AssetStorage()

        # 2. Handle 'keep_final' logic
        if keep_final:
            # Check if there's any Annotation snapshot
            latest_annotation = AnalysisSnapshotsDAO.get_latest_snapshot(dataset_id, "Annotation")

            if latest_annotation:
                # Keep the latest Annotation, delete all other snapshots
                all_snapshots = AnalysisSnapshotsDAO.get_snapshots_by_dataset(dataset_id)
                deleted_count = 0

                for snap in all_snapshots:
                    # Skip deletion for the latest annotation snapshot
                    if snap.snapshot_id != latest_annotation.snapshot_id:
                        # Determine snapshot paths
                        snap_dir_key = f"projects/{project_id}/{dataset_id}/{snap.branch_name}/{snap.snapshot_id}"
                        snap_file_path = str(snap.snapshot_path)

                        # Delete from File System
                        storage.delete_path(snap_dir_key)
                        storage.delete_path(snap_file_path)

                        # Delete from Database
                        AnalysisSnapshotsDAO.delete_snapshot(snap.id)
                        deleted_count += 1

                return {
                    "msg": f"Dataset {dataset_id} cleaned. Kept Annotation snapshot '{latest_annotation.snapshot_id}'. Deleted {deleted_count} obsolete snapshots."
                }

            # If keep_final is True but no Annotation snapshot exists, log and fall through to delete the entire dataset
            logger.info(f"No Annotation snapshot found for dataset {dataset_id}. Proceeding to delete entire dataset.")

        # 3. Full Deletion: Delete entire dataset from Database
        db_success = DatasetDAO.delete_dataset(dataset.id)

        if not db_success:
            raise RuntimeError("Failed to delete dataset record from database.")

        # 4. Full Deletion: Delete entire dataset from File System
        io_success = storage.delete_path(dataset_dir_key)

        if not io_success:
            logger.warning(f"Dataset record deleted, but failed to cleanup files at: {dataset_dir_key}")

        return {"msg": f"Dataset {dataset_id} and all associated files deleted successfully"}
