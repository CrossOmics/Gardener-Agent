from typing import Optional

from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from infrastructure.database.dao.dataset_dao import DatasetDAO


class ProjectIdResolverService:
    def __init__(self):
        pass

    def resolve_project_id(
            self,
            project_id: Optional[str] = None,
            dataset_id: Optional[str] = None,
            snapshot_id: Optional[str] = None
    ) -> str:
        """
        Ensures a valid project_id exists.
        If project_id is None, it attempts to find it via dataset_id or snapshot_id.
        """

        # If project_id is already provided, return it immediately.
        if project_id is not None and project_id != 'None' and project_id.startswith("p"):
            return project_id

        # Try to resolve via dataset_id
        if dataset_id:
            dataset = DatasetDAO.get_dataset_by_business_id(dataset_id)
            if dataset:
                # First project_id is an Object
                return dataset.project_id.project_id

        # Try to resolve via snapshot_id
        if snapshot_id:
            # Use the DAO method we defined earlier
            snapshot = AnalysisSnapshotsDAO.get_snapshot_by_business_id(snapshot_id)

            if snapshot and snapshot.dataset_id:
                # First project_id is an Object
                return snapshot.dataset_id.project_id.project_id

        # if still can't find it, just return None
        return None

    def resolve_dataset_id(self, snapshot_id: str) -> Optional[str]:
        """
        Resolves the dataset_id associated with a given snapshot_id.
        Useful when the request only provides the snapshot ID.
        """
        if not snapshot_id:
            return None

        snapshot = AnalysisSnapshotsDAO.get_snapshot_by_business_id(snapshot_id)

        # Check if snapshot exists and has a linked dataset
        if snapshot and snapshot.dataset_id:
            # Return the business ID of the dataset
            return snapshot.dataset_id.dataset_id

        return None
