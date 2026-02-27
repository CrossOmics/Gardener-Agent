from typing import List
from fastapi import APIRouter, Depends, HTTPException, status

from core.task_executor import task_executor
from service.search_annotation_methods_service import AnnotationMethodService
from dto.request.annotation_method_request import SearchAnnotationMethodRequest
from dto.response.annotation_method_response import AnnotationMethodResponse

router = APIRouter(
    prefix="/api/v1/methods/annotation",
    tags=["Annotation Reference Methods"]
)


@router.get(
    "/search",
    response_model=List[AnnotationMethodResponse],
    summary="Search Annotation Methods",
    description="Search for reference methods (GSEApy libraries or CellTypist models) by keyword and type."
)
async def search_annotation_methods(
        request: SearchAnnotationMethodRequest = Depends(),
        service: AnnotationMethodService = Depends()
):
    """
    Search endpoint.

    Args:
        request: Query parameters automatically mapped to SearchAnnotationMethodRequest DTO.
        service: Injected AnnotationMethodService.

    Returns:
        List[AnnotationMethodResponse]: A list of matching annotation methods.
    """
    def _task():
        return service.search_methods(request)

    try:
        return await task_executor.run_in_thread(_task)
    except Exception as e:
        # Log the error here in a real application
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"An error occurred while searching annotation methods: {str(e)}"
        )
