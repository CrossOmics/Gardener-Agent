from typing import Dict

from pydantic import BaseModel, Field


class DisplayNameResponse(BaseModel):
    """
    Data Transfer Object returning a map of IDs to their display names.
    """
    id_type: str = Field(..., description="The category of the entities.")
    names: Dict[str, str] = Field(
        ...,
        description="A dictionary mapping the input current_id to its resolved display name."
    )
