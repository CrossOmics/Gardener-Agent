
from fastapi import APIRouter, Depends, HTTPException, status
from loguru import logger

from core.task_executor import task_executor
from service.annotation_service import AnnotationService
from dto.request.annotation_request import RunAnnotationRequest, RunCellTypistRequest, RunFullAnnotationRequest, \
    UpdateAnnotationLabelRequest
from dto.response.annotation_result_dto import AnnotationResultDTO, UpdateAnnotationLabelResponse

router = APIRouter(
    prefix="/api/v1/annotation",
    tags=["Biological Annotation"]
)


@router.post(
    "/full",
    response_model=AnnotationResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Run Full Annotation Pipeline",
    description="Executes a consolidated annotation pipeline. Can run GSEApy (Enrichr), CellTypist, or both, based on the provided parameters. Saves all results into a single snapshot."
)
async def run_full_annotation_pipeline(
        request: RunFullAnnotationRequest,
        annotation_service: AnnotationService = Depends()
):
    """
    Executes the unified annotation pipeline.
    """
    try:
        logger.info(f"Received Full Annotation request for Project {request.project_id}")
        return await task_executor.run_in_thread(annotation_service.run_full_annotation, request)
    except ValueError as e:
        logger.warning(f"Annotation validation error: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        logger.error(f"Internal error during full annotation: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Annotation pipeline failed: {str(e)}")


@router.patch(
    "/labels/update",
    response_model=UpdateAnnotationLabelResponse,
    status_code=status.HTTP_200_OK,
    summary="Update Annotation Labels",
    description="Manually updates the predicted cell type labels for specific clusters in an existing annotation snapshot."
)
async def update_annotation_labels(
        request: UpdateAnnotationLabelRequest,
        annotation_service: AnnotationService = Depends()
):
    """
    Updates cluster labels in the snapshot metadata.
    """
    try:
        logger.info(f"Received label update request for Snapshot {request.snapshot_id}")
        return await task_executor.run_in_thread(annotation_service.update_cluster_labels, request)
    except ValueError as e:
        error_msg = str(e)
        if "not found" in error_msg.lower():
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_msg)
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)
    except Exception as e:
        logger.error(f"Internal error during label update: {e}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"Failed to update labels: {str(e)}")


@router.post("/gseapy/all", response_model=AnnotationResultDTO, status_code=status.HTTP_201_CREATED)
async def run_all_gseapy_cluster_annotation(
        request: RunAnnotationRequest,
        annotation_service: AnnotationService = Depends()
):
    """
    Run automated batch annotation for ALL clusters using GSEApy.
    """
    try:
        logger.info(f"Received BATCH annotation request for Project {request.project_id}")
        return await task_executor.run_in_thread(annotation_service.run_all_gseapy_annotation, request)
    except ValueError as e:
        logger.warning(f"Batch annotation validation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Internal error during batch annotation: {e}")
        raise HTTPException(status_code=500, detail=f"Batch annotation failed: {str(e)}")


@router.post("/celltypist", response_model=AnnotationResultDTO, status_code=status.HTTP_201_CREATED)
async def run_celltypist_auto_annotation(
        request: RunCellTypistRequest,
        annotation_service: AnnotationService = Depends()
):
    """
    Run automated cell type annotation (CellTypist) for the entire dataset.
    """
    try:
        logger.info(f"Received CellTypist request for Project {request.project_id} using model {request.model_names}")
        return await task_executor.run_in_thread(annotation_service.run_celltypist_annotations, request)
    except ValueError as e:
        logger.warning(f"CellTypist validation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Internal error during CellTypist annotation: {e}")
        raise HTTPException(status_code=500, detail=f"CellTypist annotation failed: {str(e)}")


@router.post("/gseapy/single", response_model=AnnotationResultDTO, status_code=status.HTTP_201_CREATED)
async def run_single_gseapy_cluster_annotation(
        request: RunAnnotationRequest,
        annotation_service: AnnotationService = Depends()
):
    """
    Run functional enrichment analysis (GSEApy/Enrichr) for a SINGLE cluster.
    """
    try:
        logger.info(f"Received annotation request for Cluster {request.cluster_id} in Project {request.project_id}")
        return await task_executor.run_in_thread(annotation_service.run_single_gseapy_annotation, request)
    except ValueError as e:
        logger.warning(f"Annotation validation failed: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.error(f"Internal error during annotation: {e}")
        raise HTTPException(status_code=500, detail=f"Annotation failed: {str(e)}")
