from typing import List

from fastapi import APIRouter, Depends, status, Path, HTTPException, Query

from dto.response.snapshot_ancestor_node import SnapshotAncestorNode
from dto.response.snapshot_query_response import SnapshotQueryResponse
from service.snapshot_service import SnapshotService
from dto.request.snapshot_dto import CreateSnapshotRequest, UpdateSnapshotRequest
from dto.response.snapshot_response import SnapshotResponse

router = APIRouter(
    prefix="/api/v1/snapshots",
    tags=["Analysis Snapshots"]
)


@router.get(
    "/query/{snapshot_id}",
    response_model=SnapshotQueryResponse,
    summary="Query Snapshot Details (Parameters & Metrics)",
    description="Retrieves experiment details. Returns 'input' (parameters) and 'output' (biological metrics), automatically excluding heavy file paths."
)
def query_snapshot_detail(
        snapshot_id: str = Path(..., description="The unique business ID of the snapshot"),
        snapshot_service: SnapshotService = Depends()
):
    """
    Endpoint to fetch specific snapshot details.
    Used by the frontend to display the 'Parameters' and 'Results' tabs for a specific experiment node.
    """
    try:
        return snapshot_service.query_snapshot(snapshot_id)

    except ValueError as e:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=str(e)
        )
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to query snapshot: {str(e)}"
        )


@router.get(
    "/stage/query/",
    response_model=SnapshotQueryResponse,
    summary="Query Latest Snapshot by Stage",
    description="Retrieves the most recent snapshot for a specific dataset and analysis branch (e.g., 'Preprocessing'). Returns inputs and outputs."
)
def query_latest_stage_snapshot(
        dataset_id: str = Query(..., description="The business ID of the dataset"),
        branch_name: str = Query(..., description="The branch name (e.g., 'Preprocessing', 'Clustering')"),
        service: SnapshotService = Depends()
):
    """
    Endpoint to fetch the latest snapshot for a specific analysis stage.
    Useful for retrieving the current state of a step without knowing the specific Snapshot ID.
    """
    try:
        return service.query_latest_snapshot(dataset_id, branch_name)
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to query latest snapshot: {str(e)}"
        )


@router.post("/create", response_model=SnapshotResponse, status_code=status.HTTP_201_CREATED)
def create_snapshot(
        request: CreateSnapshotRequest,
        service: SnapshotService = Depends()
):
    """
    Create a new analysis snapshot manually.
    """
    return service.create_snapshot(request)


@router.patch("/update/{snapshot_id}", response_model=SnapshotResponse)
def update_snapshot(
        snapshot_id: str,
        request: UpdateSnapshotRequest,
        service: SnapshotService = Depends()
):
    """
    Update snapshot metadata (e.g., add user notes or rename branch).
    """
    return service.update_snapshot(snapshot_id, request)


@router.delete("/delete/{snapshot_id}", status_code=status.HTTP_200_OK)
def delete_snapshot(
        snapshot_id: str,
        service: SnapshotService = Depends()
):
    """
    Delete a snapshot record.
    """
    return service.delete_snapshot(snapshot_id)


@router.delete(
    "/stage/delete",
    status_code=status.HTTP_200_OK,
    summary="Delete All Snapshots in a Stage",
    description="Deletes all snapshots associated with a specific dataset and analysis stage (branch_name). Optionally keeps the most recent one."
)
def delete_snapshots_by_stage(
        dataset_id: str = Query(..., description="The business ID of the dataset"),
        stage_name: str = Query(..., description="The analysis stage/branch name (e.g., 'Preprocessing')"),
        keep_latest: bool = Query(False,
                                  description="If true, the most recent snapshot in this stage will not be deleted"),
        service: SnapshotService = Depends()
):
    """
    Endpoint to clear out snapshots for a specific pipeline step.
    Useful for resetting a workflow stage or clearing clutter.
    """
    try:
        result = service.delete_snapshots_by_stage(
            dataset_id=dataset_id,
            stage_name=stage_name,
            keep_latest=keep_latest
        )
        return result
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to delete snapshots for stage '{stage_name}': {str(e)}"
        )


@router.get(
    "/lineage/ancestors/{snapshot_id}",
    response_model=List[SnapshotAncestorNode],
    status_code=status.HTTP_200_OK,
    summary="Get Snapshot Ancestors",
    description="Recursively traces the parent lineage of a given snapshot up to the root."
)
def query_snapshot_ancestors(
        snapshot_id: str = Path(..., description="The unique business ID of the target snapshot"),
        service: SnapshotService = Depends()
):
    """
    Endpoint to fetch the upward lineage path of a specific snapshot.
    """
    try:
        return service.get_snapshot_ancestors(snapshot_id)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to trace lineage for snapshot '{snapshot_id}': {str(e)}"
        )


@router.get(
    "/lineage/ancestors/full/{snapshot_id}",
    response_model=List[SnapshotQueryResponse],
    status_code=status.HTTP_200_OK,
    summary="Get Full Snapshot Ancestors",
    description="Recursively traces the parent lineage of a given snapshot up to the root, returning full details for each node."
)
def query_snapshot_ancestors_full(
        snapshot_id: str = Path(..., description="The unique business ID of the target snapshot"),
        service: SnapshotService = Depends()
):
    """
    Endpoint to fetch the upward lineage path of a specific snapshot with full details.
    """
    try:
        return service.get_snapshot_ancestors_full(snapshot_id)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to trace lineage for snapshot '{snapshot_id}': {str(e)}"
        )


@router.get(
    "/peers/{snapshot_id}",
    response_model=List[SnapshotQueryResponse],
    status_code=status.HTTP_200_OK,
    summary="Get Peer Snapshots",
    description="Retrieves all other snapshots that share the same dataset and branch as the given snapshot."
)
def query_snapshot_peers(
        snapshot_id: str = Path(..., description="The unique business ID of the reference snapshot"),
        service: SnapshotService = Depends()
):
    """
    Endpoint to fetch sibling/peer snapshots in the same branch.
    """
    try:
        return service.get_snapshot_peers(snapshot_id)
    except HTTPException as e:
        raise e
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to query peer snapshots: {str(e)}"
        )
