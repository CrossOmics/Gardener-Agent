
from typing import Annotated
from fastapi import APIRouter, Depends, status, HTTPException
from loguru import logger

from core.task_executor import task_executor
from dto.request.merge_clusters_request import MergeClustersRequest
from dto.request.run_clustering_request import RunClusteringRequest
from dto.request.sub_clustering_request import SubClusteringRequest
from dto.response.clustering_result_dto import ClusteringResultDTO
from service.clustering_service import ClusteringService

# Define Router
router = APIRouter(prefix="/api/v1/clustering", tags=["Analysis"])

# Dependency Injection
ServiceDep = Annotated[ClusteringService, Depends()]


@router.post(
    "/create",
    response_model=ClusteringResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Run Clustering & UMAP",
    description="Performs community detection (Leiden/Louvain) at multiple resolutions and computes UMAP embeddings.",
    responses={
        400: {"description": "Invalid resolution settings or unsupported method"},
        404: {"description": "Source snapshot not found"},
        500: {"description": "Clustering calculation failed (e.g., missing neighbors)"}
    }
)
async def run_clustering(
        request: RunClusteringRequest,
        service: ServiceDep
):
    """
    Clustering & UMAP Projection.
    """
    try:
        logger.info(f"Starting clustering ({request.method}) for project: {request.project_id}")
        return await task_executor.run_in_thread(service.run_clustering, request)

    except ValueError as e:
        error_msg = str(e)
        logger.error(f"Clustering validation error: {error_msg}")

        if "not found" in error_msg.lower():
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_msg)
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)

    except Exception as e:
        # Handle unexpected runtime errors
        logger.exception("Unexpected error during clustering execution")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CLUSTERING_ERROR: {str(e)}"
        )


@router.post(
    "/merge",
    response_model=ClusteringResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Merge Clusters",
    description="Merges specified clusters into a single new group, updates the AnnData object, and regenerates visualizations (UMAP/Dendrogram) in a new snapshot.",
    responses={
        400: {"description": "Invalid cluster IDs or parameters"},
        404: {"description": "Source snapshot not found"},
        500: {"description": "Merge operation failed"}
    }
)
async def merge_clusters(
        request: MergeClustersRequest,
        service: ServiceDep
):
    """
    Merge specific clusters manually.
    """
    try:
        logger.info(f"Merging clusters for project: {request.project_id}, snapshot: {request.snapshot_id}")
        return await task_executor.run_in_thread(service.merge_clusters, request)

    except ValueError as e:
        # Handle validation errors (e.g., cluster ID not found, snapshot missing)
        error_msg = str(e)
        logger.error(f"Cluster merge validation error: {error_msg}")

        if "not found" in error_msg.lower():
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_msg)
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)

    except Exception as e:
        # Handle unexpected system errors
        logger.exception("Unexpected error during cluster merge")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"MERGE_ERROR: {str(e)}"
        )


@router.post(
    "/sub",
    response_model=ClusteringResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Run Sub-Clustering",
    description="Performs hierarchical sub-clustering on specific target clusters. Can overwrite the existing snapshot or create a new branch.",
    responses={
        400: {"description": "Invalid target clusters or method parameters"},
        404: {"description": "Source snapshot or target cluster not found"},
        500: {"description": "Sub-clustering calculation failed"}
    }
)
async def run_sub_clustering(
        request: SubClusteringRequest,
        service: ServiceDep
):
    """
    Sub-Cluster specific groups.
    """
    try:
        logger.info(f"Sub-clustering targets {request.target_clusters} for snapshot: {request.snapshot_id}")
        return await task_executor.run_in_thread(service.run_sub_clustering, request)

    except ValueError as e:
        # Handle validation errors (e.g., target cluster not found)
        error_msg = str(e)
        logger.error(f"Sub-clustering validation error: {error_msg}")

        if "not found" in error_msg.lower():
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_msg)
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)

    except Exception as e:
        # Handle unexpected system errors
        logger.exception("Unexpected error during sub-clustering")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"SUB_CLUSTERING_ERROR: {str(e)}"
        )
