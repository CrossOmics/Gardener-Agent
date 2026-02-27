
from fastapi import APIRouter, Depends, HTTPException, status, Query
from loguru import logger

from core.task_executor import task_executor
from service.user_preference_service import UserPreferenceService
from dto.request.user_preference_request import CreateUserPreferenceRequest, UpdateUserPreferenceRequest
from dto.response.user_preference_response import UserPreferenceResponse

router = APIRouter(
    prefix="/api/v1/user/preference",
    tags=["User Preferences"]
)


@router.post(
    "/create",
    response_model=UserPreferenceResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Create New User Preference",
    description="Creates a new global user preference setting. If no name is provided, a default one is generated."
)
async def create_user_preference(
        request: CreateUserPreferenceRequest,
        service: UserPreferenceService = Depends()
):
    """
    Creates a new user preference entry.

    Workflow:
    1. Validates the input request.
    2. Generates a default name if 'preference_name' is missing.
    3. Persists the settings JSON to the database.
    4. Returns the created preference object with timestamps.
    """
    def _task():
        logger.info(f"Creating new user preference: {request.preference_name or 'Auto-generated'}")
        return service.create_preference(request)

    try:
        return await task_executor.run_in_thread(_task)
    except HTTPException as e:
        # Re-raise HTTP exceptions from service (e.g., duplicate name)
        raise e
    except Exception as e:
        logger.error(f"Internal error during preference creation: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to create preference: {str(e)}"
        )


@router.get(
    "/query",
    response_model=UserPreferenceResponse,
    summary="Get User Preference by Name",
    description="Retrieves a specific user preference by its unique name. Also updates the 'updated_at' timestamp to mark it as recently accessed."
)
async def get_user_preference(
        name: str = Query(..., description="Unique name of the preference"),
        service: UserPreferenceService = Depends()
):
    """
    Retrieves a preference and touches its timestamp.

    Workflow:
    1. Searches for the preference by name (ignoring soft-deleted ones).
    2. If found, updates the 'updated_at' field to the current time.
    3. Returns the preference details.
    """
    def _task():
        return service.get_preference(name)

    try:
        return await task_executor.run_in_thread(_task)
    except HTTPException as e:
        raise e
    except Exception as e:
        logger.error(f"Internal error retrieving preference '{name}': {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve preference: {str(e)}"
        )


@router.get(
    "/latest/one",
    response_model=UserPreferenceResponse,
    summary="Get Latest User Preference",
    description="Fetches the single most recently created or updated user preference."
)
async def get_latest_preference(
        service: UserPreferenceService = Depends()
):
    """
    Fetches the most recent preference.

    Workflow:
    1. Queries the database for the record with the latest 'updated_at' timestamp.
    2. Returns the single latest preference object.
    """
    def _task():
        return service.get_latest_preference()

    try:
        return await task_executor.run_in_thread(_task)
    except HTTPException as e:
        raise e
    except Exception as e:
        logger.error(f"Internal error retrieving latest preference: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to fetch latest preference: {str(e)}"
        )


@router.patch(
    "/update",
    response_model=UserPreferenceResponse,
    summary="Update User Preference",
    description="Updates the settings JSON and optionally renames an existing preference."
)
async def update_user_preference(
        request: UpdateUserPreferenceRequest,
        service: UserPreferenceService = Depends()
):
    """
    Updates an existing preference.

    Workflow:
    1. Locates the preference by the current 'name'.
    2. Updates the 'settings' JSON string.
    3. If 'new_name' is provided and different, renames the preference.
    4. Refreshes the 'updated_at' timestamp.
    """
    def _task():
        logger.info(f"Updating preference '{request.name}'")
        return service.update_preference(request)

    try:
        return await task_executor.run_in_thread(_task)
    except HTTPException as e:
        raise e
    except Exception as e:
        logger.error(f"Internal error updating preference '{request.name}': {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to update preference: {str(e)}"
        )


@router.delete(
    "/delete",
    status_code=status.HTTP_204_NO_CONTENT,
    summary="Delete User Preference",
    description="Soft-deletes a user preference by name. The record is marked as deleted but not physically removed."
)
async def delete_user_preference(
        name: str = Query(..., description="Unique name of the preference to delete"),
        service: UserPreferenceService = Depends()
):
    """
    Soft deletes a preference.

    Workflow:
    1. Locates the preference by name.
    2. Sets the 'is_deleted' flag to True.
    3. Returns 204 No Content on success.
    """
    def _task():
        logger.info(f"Deleting user preference: {name}")
        service.delete_preference(name)
        return

    try:
        return await task_executor.run_in_thread(_task)
    except HTTPException as e:
        raise e
    except Exception as e:
        logger.error(f"Internal error deleting preference '{name}': {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to delete preference: {str(e)}"
        )
