from typing import Optional, Dict, Any
from pydantic import BaseModel, Field


class CreateUserPreferenceRequest(BaseModel):
    """
    Request payload for creating a new user preference.
    """
    preference_name: Optional[str] = Field(
        None,
        description="Optional custom name. If empty, a default name is generated."
    )
    settings: Dict[str, Any] = Field(..., description="Pipeline parameters object.")


class UpdateUserPreferenceRequest(BaseModel):
    """
    Request payload for updating an existing user preference.
    Contains the target name, new settings, and an optional new name.
    """
    name: str = Field(
        ...,
        description="The current unique name of the preference to update."
    )
    settings: Dict[str, Any] = Field(..., description="The new pipeline parameters object.")
    new_name: Optional[str] = Field(
        None,
        description="Optional new name. If provided, the preference will be renamed."
    )
