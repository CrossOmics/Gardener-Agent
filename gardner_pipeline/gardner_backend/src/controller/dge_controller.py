
from typing import Annotated
from fastapi import APIRouter, Depends, status, HTTPException
from loguru import logger

from core.task_executor import task_executor
from dto.request.run_dge_request import RunDGERequest
from dto.response.dge_result_dto import DGEResultDTO
from service.dge_service import DGEService

router = APIRouter(prefix="/api/v1/dge", tags=["DGE"])
ServiceDep = Annotated[DGEService, Depends()]


@router.post(
    "/new",
    response_model=DGEResultDTO,
    status_code=status.HTTP_201_CREATED,
    summary="Differential Expression Gene (DGE) Analysis",
    description="Identifies marker genes for each cluster using statistical tests (Wilcoxon/t-test) and generates visualization plots.",
    responses={
        400: {"description": "Invalid groupby column or insufficient groups (<2)"},
        404: {"description": "Source snapshot not found"},
        500: {"description": "Calculation failed (e.g., singular matrix, data integrity issues)"}
    }
)
async def run_deg_analysis(
        request: RunDGERequest,
        service: ServiceDep
):
    """
    Differential Expression Analysis.
    """
    try:
        logger.info(f"Starting DGE analysis for project: {request.project_id} on group: {request.groupby}")
        return await task_executor.run_in_thread(service.run_dge_analysis, request)

    except ValueError as e:
        # Handle logical errors (e.g., column missing, single group)
        error_msg = str(e)
        logger.error(f"DGE validation error: {error_msg}")

        if "not found" in error_msg.lower() and "column" not in error_msg.lower():
            # Likely a snapshot not found error
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=error_msg)
        else:
            # Likely a groupby column missing or insufficient groups error
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=error_msg)

    except Exception as e:
        # Handle unexpected runtime errors
        logger.exception("Unexpected error during DGE execution")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"CALCULATION_ERROR: {str(e)}"
        )
