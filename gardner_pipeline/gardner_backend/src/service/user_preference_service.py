from typing import Optional
from fastapi import Depends, HTTPException, status

from infrastructure.database.dao.user_preference_dao import UserPreferenceDAO
from dto.request.user_preference_request import CreateUserPreferenceRequest, UpdateUserPreferenceRequest
from dto.response.user_preference_response import UserPreferenceResponse
from infrastructure.database.model.default_preference import standard_settings


class UserPreferenceService:
    """
    Service layer for User Preference logic.
    """

    def __init__(self, dao: UserPreferenceDAO = Depends()):
        self.dao = dao

    def get_preference(self, name: str) -> UserPreferenceResponse:
        """
        1. Retrieve preference by name.
        2. Update the 'updated_at' timestamp (touch) to mark it as recently used.
        3. Return the response DTO.
        """
        # 1. Fetch from DB
        pref = self.dao.get_by_name(name)
        if not pref:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Preference '{name}' not found."
            )
        # 2. Touch timestamp (updates the last access time)
        self.dao.update_timestamp(name)
        # 3. Convert to DTO
        return UserPreferenceResponse.model_validate(pref)

    def create_preference(self, request: CreateUserPreferenceRequest) -> UserPreferenceResponse:
        """
        Create a new preference.
        Logic:
        - If preference_name is empty, generate "preference_{count + 1}".
        - Persist to DB.
        - Return DTO.
        """
        name = request.preference_name.strip() if request.preference_name else ""
        # 1. Handle auto-naming
        if not name:
            count = self.dao.count_all()
            name = f"preference_{count + 1}"
        # 2. Create in DB
        new_pref = self.dao.create_preference(name, request.settings)
        if not new_pref:
            # Creation failed, likely due to duplicate name
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Failed to create preference. Name '{name}' might already exist."
            )
        # 3. Convert to DTO
        return UserPreferenceResponse.model_validate(new_pref)

    def get_latest_preference(self) -> Optional[UserPreferenceResponse]:
        """
        Fetch the most recently updated/created preference.
        """
        pref = self.dao.get_latest_updated()
        if not pref:
            count = self.dao.count_all()
            # Fallback: No preference found, initialize a default one
            default_name = f"preference_{count + 1}"
            # Create and persist the default preference
            pref = self.dao.create_preference(default_name, standard_settings)
            # Validation check
            if not pref:
                raise HTTPException(
                    status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                    detail="Failed to create default preference."
                )
        return UserPreferenceResponse.model_validate(pref)

    def delete_preference(self, name: str) -> None:
        """
        Soft delete a preference by name.
        """
        success = self.dao.delete_by_name(name)
        if not success:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Preference '{name}' not found or already deleted."
            )

    def update_preference(self, request: UpdateUserPreferenceRequest) -> UserPreferenceResponse:
        """
        Update settings and optionally rename an existing preference.
        Uses ID-based retrieval for robustness.
        """
        # 1. Attempt update and get the ID back
        updated_id = self.dao.update_preference(request.name, request.settings, request.new_name)
        if updated_id == -1:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail=f"Update failed. Preference '{request.name}' not found or new name '{request.new_name}' already exists."
            )
        # 2. Retrieve fresh object directly by ID
        # This avoids any confusion about which name to query (old vs new)
        updated_pref = self.dao.get_by_id(updated_id)
        if not updated_pref:
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail="Error retrieving updated preference data."
            )
        # 3. Return DTO
        return UserPreferenceResponse.model_validate(updated_pref)
