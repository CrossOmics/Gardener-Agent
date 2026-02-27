from typing import List, Optional, Dict, Any
from fastapi import Depends, HTTPException

from dto.response.snapshot_query_response import SnapshotQueryResponse
from infrastructure.database.dao.analysis_snapshots_dao import AnalysisSnapshotsDAO
from dto.request.snapshot_dto import CreateSnapshotRequest, UpdateSnapshotRequest
from infrastructure.filesystem.storage import AssetStorage
from service.helper.project_id_resolver_service import ProjectIdResolverService
from loguru import logger


class SnapshotService:
    def __init__(self, dao: AnalysisSnapshotsDAO = Depends(),
                 storage: AssetStorage = Depends(),
                 project_id_resolver: ProjectIdResolverService = Depends()):
        self.dao = dao
        self.storage = storage
        self.project_id_resolver = project_id_resolver

    def _construct_snapshot_response(self, snapshot) -> SnapshotQueryResponse:
        """
        Private helper to format the snapshot object into the API response format.
        Handles the logic of separating params (input) and filtering metrics (output).
        """
        # 1. Prepare 'input' (Experiment Parameters)
        input_data = snapshot.params_json if snapshot.params_json else {}

        # 2. Prepare 'output' (Metrics only)
        # Filter thumbnail_json to exclude file paths and plot references
        raw_output = snapshot.thumbnail_json if snapshot.thumbnail_json else {}
        filtered_output = {}

        # Keywords to identify file assets to exclude
        exclusion_keywords = ["plot", "path", "json", "csv"]

        for key, value in raw_output.items():
            # Convert key to lowercase for case-insensitive matching
            key_lower = key.lower()

            # Keep the key ONLY if it does not contain any exclusion keyword
            if not any(keyword in key_lower for keyword in exclusion_keywords):
                filtered_output[key] = value

        # 3. Assemble and return the DTO
        return SnapshotQueryResponse(
            snapshot_id=snapshot.snapshot_id,
            parent_snapshot_id=snapshot.parent_snapshot_id,
            branch_name=snapshot.branch_name,
            create_time=snapshot.create_time,
            snapshot_name=snapshot.snapshot_name,
            user_notes=snapshot.user_notes,
            input=input_data,
            output=filtered_output
        )

    def get_snapshot(self, snapshot_id: str):
        """
        Fetch a single snapshot by its business ID (string UUID).
        """
        snapshot = self.dao.get_snapshot_by_business_id(snapshot_id)
        if not snapshot:
            raise HTTPException(status_code=404, detail=f"Snapshot {snapshot_id} not found")
        return snapshot

    def query_snapshot(self, snapshot_id: str) -> SnapshotQueryResponse:
        """
        Retrieves snapshot details by specific Snapshot ID.
        """
        # Retrieve the snapshot object (raises 404 if not found)
        snapshot = self.get_snapshot(snapshot_id)

        # Use helper to format response
        return self._construct_snapshot_response(snapshot)

    def query_latest_snapshot(self, dataset_id: str, branch_name: str) -> SnapshotQueryResponse:
        """
        Retrieves the LATEST snapshot details for a specific Dataset and Branch.
        Useful for fetching the most recent 'Preprocessing' or 'Clustering' result without knowing the ID.
        """
        # 1. Query DAO for the latest record
        snapshot = self.dao.get_latest_snapshot(dataset_id, branch_name)

        if not snapshot:
            raise HTTPException(
                status_code=404,
                detail=f"No snapshot found for dataset {dataset_id} on branch '{branch_name}'"
            )

        # 2. Use helper to format response
        return self._construct_snapshot_response(snapshot)

    def list_snapshots_by_dataset(self, dataset_id: str):
        """
        List all snapshots belonging to a specific dataset.
        """
        return self.dao.get_snapshots_by_dataset(dataset_id)

    def create_snapshot(self, request: CreateSnapshotRequest):
        """
        Create a new snapshot record.
        """
        # Check if ID already exists to avoid conflict
        existing = self.dao.get_snapshot_by_business_id(request.snapshot_id)
        if existing:
            raise HTTPException(status_code=400, detail="Snapshot ID already exists")

        snapshot = self.dao.create_snapshot(
            dataset_id=request.dataset_id,
            snapshot_id=request.snapshot_id,
            branch_name=request.branch_name,
            snapshot_path=request.snapshot_path,
            parent_snapshot_id=request.parent_snapshot_id,
            params_json=request.params_json,
            thumbnail_json=request.thumbnail_json,
            user_notes=request.user_notes
        )

        if not snapshot:
            raise HTTPException(status_code=500, detail="Failed to create snapshot")

        return snapshot

    def update_snapshot(self, snapshot_id: str, request: UpdateSnapshotRequest):
        """
        Update metadata (notes, branch name) for an existing snapshot.
        """
        # 1. Find the snapshot first to get its Primary Key (ID)
        snapshot = self.get_snapshot(snapshot_id)

        # 2. Call DAO update using the integer PK
        success = self.dao.update_snapshot(
            pk_id=snapshot.id,
            branch_name=request.branch_name,
            user_notes=request.user_notes,
            params_update=request.params_update
        )

        if not success:
            raise HTTPException(status_code=500, detail="Failed to update snapshot")

        return self.dao.get_snapshot_by_business_id(snapshot.id)

    def delete_snapshot(self, snapshot_id: str):
        """
        Deletes a snapshot record from the database AND its corresponding files from storage.
        """
        # 1. Retrieve the snapshot object (ensures existence and gets the integer PK)
        snapshot = self.dao.get_snapshot_by_business_id(snapshot_id)
        if not snapshot:
            raise HTTPException(status_code=404, detail=f"Snapshot {snapshot_id} not found")

        # 2. Resolve IDs to construct the file path
        project_id = self.project_id_resolver.resolve_project_id(snapshot_id=snapshot_id)
        dataset_id = self.project_id_resolver.resolve_dataset_id(snapshot_id=snapshot_id)

        # 3. Determine the relative directory path for this snapshot

        # Get the parent directory of the .h5ad file
        snapshot_path = str(snapshot.snapshot_path)
        # Fallback construction (less reliable if paths vary)
        relative_dir_key = f"projects/{project_id}/{dataset_id}/{snapshot.branch_name}/{snapshot_id}"

        # 4. Delete from Database
        db_success = self.dao.delete_snapshot(snapshot.id)
        if not db_success:
            raise HTTPException(status_code=500, detail="Failed to delete snapshot record from database")

        # 5. Delete from File System
        io_success = self.storage.delete_path(relative_dir_key) and self.storage.delete_path(snapshot_path)

        if not io_success:
            logger.warning(f"Snapshot record deleted, but failed to cleanup files at: {relative_dir_key}")

        return {"msg": f"Snapshot {snapshot_id} deleted successfully"}

    def get_lineage(self, dataset_id: str) -> List[Dict]:
        """
        Get the parent-child relationship tree for the dataset.
        """
        return self.dao.get_snapshot_lineage_by_dataset(dataset_id)

    def update_snapshot_name(self, snapshot_id: str, new_name: str) -> bool:
        """
        Updates the user-friendly display name of a specific snapshot.
        Errors are handled at the controller layer.
        """
        # Verify snapshot exists first
        self.get_snapshot(snapshot_id)

        # 2. Call the DAO method to perform the database update
        success = self.dao.update_snapshot_name(snapshot_id=snapshot_id, new_name=new_name)

        # 3. Log failure if the database update did not affect any rows
        if not success:
            logger.error(f"Failed to update display name for snapshot {snapshot_id}")

        return success

    def get_snapshot_name(self, snapshot_id: str) -> str:
        """
        Retrieves the user-friendly display name of a specific snapshot.
        """
        # self.get_snapshot already raises a 404 HTTPException if not found
        snapshot = self.get_snapshot(snapshot_id)

        return snapshot.snapshot_name

    def delete_snapshots_by_stage(self, dataset_id: str, stage_name: str, keep_latest: bool = False) -> Dict[str, Any]:
        """
        Deletes all snapshots for a given dataset under a specific branch/stage.
        Optionally retains the latest snapshot.
        """
        # 1. Fetch all snapshots matching the dataset and stage
        # We assume the DAO returns them ordered by create_time descending (newest first)
        # If your DAO doesn't have a direct method for this, use the existing list_snapshots_by_dataset and filter.
        all_snapshots = self.dao.get_snapshots_by_dataset(dataset_id)

        # Filter for the specific stage
        target_snapshots = [s for s in all_snapshots if s.branch_name == stage_name]

        if not target_snapshots:
            raise ValueError(f"No snapshots found for dataset {dataset_id} on stage '{stage_name}'")

        # Sort by creation time descending to ensure the latest is at index 0
        target_snapshots.sort(key=lambda x: x.create_time, reverse=True)

        snapshots_to_delete = target_snapshots

        # 2. Handle the 'keep_latest' logic
        if keep_latest:
            if len(target_snapshots) <= 1:
                return {
                    "msg": f"Only 1 snapshot exists for stage '{stage_name}'. Kept as requested.",
                    "deleted_count": 0
                }
            # Remove the first element (the latest) from the deletion list
            snapshots_to_delete = target_snapshots[1:]

        # 3. Execute Deletion
        deleted_count = 0
        failed_count = 0

        for snapshot in snapshots_to_delete:
            try:
                # Re-use the existing atomic delete method for each snapshot
                # This ensures both DB and file system are cleaned up
                self.delete_snapshot(snapshot.snapshot_id)
                deleted_count += 1
            except Exception as e:
                logger.error(f"Failed to delete snapshot {snapshot.snapshot_id} during bulk stage deletion: {e}")
                failed_count += 1

        # 4. Formulate Response
        if failed_count > 0:
            logger.warning(f"Stage deletion partially failed. Deleted: {deleted_count}, Failed: {failed_count}")

        status_msg = f"Successfully deleted {deleted_count} snapshots for stage '{stage_name}'."
        if keep_latest:
            status_msg += " Kept the latest one."

        return {
            "msg": status_msg,
            "deleted_count": deleted_count
        }

    def get_snapshot_ancestors(self, snapshot_id: str) -> List[Dict[str, Any]]:
        """
        Traces the lineage of a snapshot upwards to the root.
        Returns an ordered list from the target node up to the root node.
        """
        lineage = []
        current_id = snapshot_id

        # Loop until we hit a root node (parent_snapshot_id is None or empty string)
        # Using a safeguard limit (e.g., 50) to prevent infinite loops in corrupted DBs
        max_depth = 50
        depth = 0

        while current_id and depth < max_depth:
            # Re-use existing DAO method to fetch the node
            node = self.dao.get_snapshot_by_business_id(current_id)

            if not node:
                if depth == 0:
                    # If the initial requested snapshot doesn't exist
                    raise ValueError(f"Snapshot {snapshot_id} not found")
                else:
                    # Parent doesn't exist (broken link), stop tracing
                    logger.warning(f"Broken lineage link: Parent {current_id} not found.")
                    break

            # Append required fields. Note: branch_name maps to stage_name
            lineage.append({
                "snapshot_id": node.snapshot_id,
                "snapshot_name": node.snapshot_name,
                "stage_name": node.branch_name,
                "create_time": node.create_time
            })

            # Move up to the parent
            current_id = node.parent_snapshot_id
            depth += 1

        if depth >= max_depth:
            logger.warning(f"Lineage trace hit max depth limit ({max_depth}) for snapshot {snapshot_id}")

        return lineage

    def get_snapshot_ancestors_full(self, snapshot_id: str) -> List[SnapshotQueryResponse]:
        """
        Traces the lineage of a snapshot upwards to the root.
        Returns an ordered list of FULL snapshot details from the target node up to the root node.
        """
        lineage = []
        current_id = snapshot_id
        max_depth = 50
        depth = 0

        while current_id and depth < max_depth:
            node = self.dao.get_snapshot_by_business_id(current_id)

            if not node:
                if depth == 0:
                    raise ValueError(f"Snapshot {snapshot_id} not found")
                else:
                    logger.warning(f"Broken lineage link: Parent {current_id} not found.")
                    break

            # Construct full response
            full_response = self._construct_snapshot_response(node)
            lineage.append(full_response)

            current_id = node.parent_snapshot_id
            depth += 1

        if depth >= max_depth:
            logger.warning(f"Lineage trace hit max depth limit ({max_depth}) for snapshot {snapshot_id}")

        return lineage

    def get_snapshot_peers(self, snapshot_id: str) -> List[SnapshotQueryResponse]:
        """
        Retrieves all snapshots that share the same dataset_id and branch_name as the given snapshot_id.
        Excludes the requested snapshot_id itself from the list.
        """
        # 1. Get the reference snapshot to find its dataset and branch
        ref_snapshot = self.get_snapshot(snapshot_id)
        dataset_id = ref_snapshot.dataset_id
        branch_name = ref_snapshot.branch_name

        # 2. Query all snapshots in that branch
        all_snapshots = self.dao.get_snapshots_by_branch(dataset_id, branch_name)

        # 3. Filter out the reference snapshot and format the rest
        peers = []
        for snap in all_snapshots:
            if snap.snapshot_id != snapshot_id:
                peers.append(self._construct_snapshot_response(snap))

        return peers
