
from fastapi import APIRouter, HTTPException, status, Depends
from typing import Annotated

from core.task_executor import task_executor
from dto.request.calculate_qc_request import CalculateQCRequest
from dto.request.filter_qc_request import FilterQCRequest
from dto.request.full_preprocessing_request import FullPreprocessingRequest
from dto.request.run_hvg_request import RunHVGRequest
from dto.response.filter_qc_response import FilterQCResponse
from dto.response.full_preprocessing_response import FullPreprocessingResponse
from dto.response.hvg_result_dto import HVGResultDTO
from dto.response.qc_result_dto import QCResultDTO
from service.preprocessing_service import PreprocessingService

router = APIRouter(prefix="/api/v1/preprocessing", tags=["Preprocessing"])
ServiceDep = Annotated[PreprocessingService, Depends()]


@router.post("/run/full", response_model=FullPreprocessingResponse, status_code=status.HTTP_201_CREATED,
             summary="Full Preprocessing Pipeline",
             description="Runs the complete preprocessing pipeline (QC -> Filter -> HVG -> PCA -> Neighbors) in memory and saves a single snapshot."
             )
async def run_full_preprocessing(
        request: FullPreprocessingRequest,
        service: ServiceDep
):
    """
    Executes the full preprocessing pipeline in a single call.
    """
    try:
        return await task_executor.run_in_thread(service.full_preprocessing, request)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"DATASET_NOT_FOUND: {str(e)}")
    except RuntimeError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"PIPELINE_ERROR: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"SYSTEM_ERROR: {str(e)}")


@router.post("/qc/calculate", response_model=QCResultDTO, status_code=status.HTTP_200_OK)
async def perform_qc_calculation(
        request: CalculateQCRequest,
        service: ServiceDep
):
    """
    Calculate QC Metrics.
    """
    try:
        return await task_executor.run_in_thread(
            service.qc_calculation,
            project_id=request.project_id,
            dataset_id=request.dataset_id,
            organism=request.organism,
            custom_prefix=request.custom_prefixes
        )
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))


@router.post("/qc/filter", response_model=FilterQCResponse, status_code=status.HTTP_201_CREATED)
async def apply_qc_filter(request: FilterQCRequest, service: ServiceDep):
    """
    Apply QC Filtering (Step 2).
    """
    try:
        return await task_executor.run_in_thread(service.apply_filter, request)
    except ValueError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except FileNotFoundError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except RuntimeError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"System Error during filtering: {str(e)}")


@router.post("/hvg", response_model=HVGResultDTO, status_code=status.HTTP_201_CREATED,
             summary="Feature Selection & Scaling",
             description="Performs Normalization -> Log1p -> HVG Identification -> Raw Backup -> Scaling."
             )
async def apply_hvg(request: RunHVGRequest, service: ServiceDep):
    """
    Feature Selection (HVG) & Scaling.
    """
    try:
        return await task_executor.run_in_thread(service.apply_hvg, request)
    except ValueError as e:
        error_msg = str(e).lower()
        if "not found" in error_msg:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"DATASET_NOT_FOUND: {str(e)}")
        else:
            raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=f"INVALID_PARAMS: {str(e)}")
    except FileNotFoundError as e:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=f"DATASET_NOT_FOUND: Physical file missing - {str(e)}")
    except RuntimeError as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"CALCULATION_ERROR: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"SYSTEM_ERROR: {str(e)}")
