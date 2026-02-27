from datetime import datetime
from typing import Dict, Any

from pydantic import BaseModel


class UserPreferenceResponse(BaseModel):
    """
    Response DTO for UserPreference.
    Excludes the 'is_deleted' field.
    """
    id: int
    preference_name: str
    settings: Dict[str, Any]
    created_at: datetime
    updated_at: datetime

    class Config:
        # Pydantic V2 configuration to allow creation from ORM objects
        from_attributes = True