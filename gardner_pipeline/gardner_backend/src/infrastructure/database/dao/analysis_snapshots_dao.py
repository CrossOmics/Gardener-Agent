from datetime import datetime
from typing import List, Optional, Dict, Any, Union
from peewee import DoesNotExist, IntegrityError

from ..model.analysis_snapshots_model import AnalysisSnapshot
from ..model.dataset_model import Dataset


class AnalysisSnapshotsDAO:
    def __init__(self):
        pass

    """
    Data Access Object (DAO) for managing AnalysisSnapshot database operations.
    """

    @staticmethod
    def get_latest_snapshot(dataset_id: str, branch_name: str) -> Optional[AnalysisSnapshot]:
        """
        Retrieves the most recently created snapshot for a specific dataset and branch.

        Args:
            dataset_id: The business ID of the dataset.
            branch_name: The specific branch/step name to filter by (e.g., 'QC Filtered').

        Returns:
            The latest AnalysisSnapshot object, or None if no match found.
        """
        try:
            return (AnalysisSnapshot
                    .select()
                    .where(
                (AnalysisSnapshot.dataset_id == dataset_id) &
                (AnalysisSnapshot.branch_name == branch_name)
            )
                    .order_by(AnalysisSnapshot.create_time.desc())  # Newest first
                    .first())
        except Exception as e:
            print(f"[DB Error] Error fetching latest snapshot: {e}")
            return None

    @staticmethod
    def create_snapshot(
            dataset_id: Union[Dataset, str],
            snapshot_id: str,
            branch_name: str,
            snapshot_path: str,
            parent_snapshot_id: Optional[str] = None,
            params_json: Optional[Dict[str, Any]] = None,
            thumbnail_json: Optional[Dict[str, Any]] = None,
            user_notes: Optional[str] = None,
            snapshot_name: str = "New Snapshot",
            create_time: Optional[datetime] = None,
    ) -> Optional[AnalysisSnapshot]:
        """
        Creates a new AnalysisSnapshot record in the database.
        Automatically strips context IDs (project_id, dataset_id, snapshot_id) from params_json.
        """
        try:
            clean_params = params_json.copy() if params_json else {}
            # 2Define keys to remove (Redundant IDs that are already stored in columns)
            keys_to_remove = ["project_id", "dataset_id", "snapshot_id"]

            # Safely remove them (pop with None default prevents errors if key is missing)
            for key in keys_to_remove:
                clean_params.pop(key, None)
            snapshot = AnalysisSnapshot.create(
                dataset_id=dataset_id,
                snapshot_id=snapshot_id,
                branch_name=branch_name,
                snapshot_path=snapshot_path,
                parent_snapshot_id=parent_snapshot_id,
                params_json=clean_params,
                thumbnail_json=thumbnail_json or {},
                user_notes=user_notes,
                snapshot_name=snapshot_name,
                create_time=create_time or datetime.utcnow()
            )
            return snapshot
        except IntegrityError as e:
            print(f"[Error] Failed to create snapshot '{snapshot_id}': {e}")
            return None
        except Exception as e:
            print(f"[Error] Unexpected error creating snapshot: {e}")
            return None

    @staticmethod
    def get_snapshot_by_id(pk_id: int) -> Optional[AnalysisSnapshot]:
        try:
            return AnalysisSnapshot.get_by_id(pk_id)
        except DoesNotExist:
            return None

    @staticmethod
    def get_snapshot_by_business_id(snapshot_id: str) -> Optional[AnalysisSnapshot]:
        try:
            return AnalysisSnapshot.get(AnalysisSnapshot.snapshot_id == snapshot_id)
        except DoesNotExist:
            return None

    @staticmethod
    def get_snapshots_by_dataset(dataset_id: str) -> List[AnalysisSnapshot]:
        try:
            query = (AnalysisSnapshot
                     .select()
                     .where(AnalysisSnapshot.dataset_id == dataset_id)
                     .order_by(AnalysisSnapshot.create_time.desc()))
            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve snapshots for dataset {dataset_id}: {e}")
            return []

    @staticmethod
    def get_snapshots_by_branch(dataset_id: str, branch_name: str) -> List[AnalysisSnapshot]:
        """
        Retrieves all snapshots for a specific dataset and branch.
        """
        try:
            query = (AnalysisSnapshot
                     .select()
                     .where(
                (AnalysisSnapshot.dataset_id == dataset_id) &
                (AnalysisSnapshot.branch_name == branch_name)
            )
                     .order_by(AnalysisSnapshot.create_time.desc()))
            return list(query)
        except Exception as e:
            print(f"[DB Error] Error fetching snapshots by branch: {e}")
            return []

    @staticmethod
    def update_snapshot(
            pk_id: int,
            branch_name: Optional[str] = None,
            user_notes: Optional[str] = None,
            end_time: Optional[datetime] = None,
            params_update: Optional[Dict[str, Any]] = None,
            thumbnail_json: Optional[Dict[str, Any]] = None,
            snapshot_name: Optional[str] = None,
    ) -> bool:
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)

            if branch_name is not None:
                snapshot.branch_name = branch_name

            if user_notes is not None:
                snapshot.user_notes = user_notes

            if end_time is not None:
                snapshot.end_time = end_time
            else:
                snapshot.end_time = datetime.utcnow()

            if params_update is not None:
                current_params = snapshot.params_json or {}
                current_params.update(params_update)
                snapshot.params_json = current_params

            if thumbnail_json is not None:
                current_thumb = snapshot.thumbnail_json or {}
                current_thumb.update(thumbnail_json)
                snapshot.thumbnail_json = current_thumb

            if snapshot_name is not None:
                snapshot.snapshot_name = snapshot_name

            snapshot.save()
            return True
        except DoesNotExist:
            print(f"[Error] Cannot update: Snapshot ID {pk_id} does not exist.")
            return False
        except Exception as e:
            print(f"[Error] Failed to update snapshot {pk_id}: {e}")
            return False

    @staticmethod
    def delete_snapshot(pk_id: int) -> bool:
        try:
            snapshot = AnalysisSnapshot.get_by_id(pk_id)
            snapshot.delete_instance()
            return True
        except DoesNotExist:
            return False
        except Exception as e:
            print(f"[Error] Failed to delete snapshot {pk_id}: {e}")
            return False

    @staticmethod
    def delete_snapshots_by_dataset(dataset_id: str) -> int:
        try:
            query = AnalysisSnapshot.delete().where(AnalysisSnapshot.dataset_id == dataset_id)
            deleted_count = query.execute()
            return deleted_count
        except Exception as e:
            print(f"[Error] Failed to batch delete snapshots: {e}")
            return 0

    @staticmethod
    def get_snapshot_lineage_by_dataset(dataset_id: str) -> List[Dict[str, str]]:
        """
        Retrieves the structural relationship (ID and Parent ID) for all snapshots
        in a dataset.

        dataset_id: The business ID of the dataset.

        returns: A list of dictionaries, e.g.:
            [
                {'snapshot_id': 'snap_1', 'parent_snapshot_id': None, 'branch_name': 'QC'},
                {'snapshot_id': 'snap_2', 'parent_snapshot_id': 'snap_1', 'branch_name': 'HVG'}
            ]
        """
        try:
            # We use .dicts() to return a lightweight dictionary instead of Model objects
            # limiting columns to only what's needed for tree construction.
            query = (AnalysisSnapshot
                     .select(
                AnalysisSnapshot.snapshot_id,
                AnalysisSnapshot.parent_snapshot_id,
                AnalysisSnapshot.branch_name,
                AnalysisSnapshot.create_time
            )
                     .where(AnalysisSnapshot.dataset_id == dataset_id)
                     .order_by(AnalysisSnapshot.create_time.asc())  # Ordered by time for logic
                     .dicts())

            return list(query)
        except Exception as e:
            print(f"[Error] Failed to retrieve lineage for dataset {dataset_id}: {e}")
            return []

    @staticmethod
    def get_parent_id_by_snapshot(snapshot_id: str) -> Optional[str]:
        """
        Helper to get just the parent ID for a specific snapshot.
        """
        try:
            snapshot = (AnalysisSnapshot
                        .select(AnalysisSnapshot.parent_snapshot_id)
                        .where(AnalysisSnapshot.snapshot_id == snapshot_id)
                        .first())
            return snapshot.parent_snapshot_id if snapshot else None
        except Exception as e:
            print(f"[Error] Failed to get parent for {snapshot_id}: {e}")
            return None

    @classmethod
    def update_snapshot_name(cls, snapshot_id: str, new_name: str) -> bool:
        """
        Updates the user-defined display name (snapshot_name) of a specific snapshot.

        Args:
            snapshot_id: The unique system identifier of the snapshot.
            new_name: The new display name provided by the user.

        Returns:
            bool: True if the update was successful and a row was changed, False otherwise.
        """
        try:
            query = AnalysisSnapshot.update(snapshot_name=new_name).where(
                AnalysisSnapshot.snapshot_id == snapshot_id
            )
            rows_updated = query.execute()
            return rows_updated > 0
        except Exception as e:
            # Replace with your actual logger
            print(f"[Error] Failed to update snapshot name for {snapshot_id}: {e}")
            return False

    @staticmethod
    def count_snapshots_by_branch(dataset_id: str, branch_name: str) -> int:
        """
        Counts the number of existing snapshots for a specific dataset and stage (branch_name).

        Args:
            dataset_id: The business ID of the dataset.
            branch_name: The name of the analysis stage (e.g., 'Clustering').

        Returns:
            int: The total count of snapshots found. Returns 0 if none exist or an error occurs.
        """
        try:
            return (AnalysisSnapshot
                    .select()
                    .where(
                (AnalysisSnapshot.dataset_id == dataset_id) &
                (AnalysisSnapshot.branch_name == branch_name)
            )
                    .count())
        except Exception as e:
            # Replace with your logger if necessary
            print(f"[DB Error] Failed to count snapshots for stage '{branch_name}': {e}")
            return 0
