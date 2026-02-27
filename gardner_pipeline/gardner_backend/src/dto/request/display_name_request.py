from pydantic import BaseModel, Field
from typing import Literal, List


class DisplayNameRequest(BaseModel):
    """
    Data Transfer Object for bulk requesting entity names.
    All IDs in 'current_ids' must belong to the specified 'id_type'.
    """
    id_type: Literal['project', 'dataset', 'snapshot'] = Field(
        ...,
        description="The category of the entities being queried."
    )
    current_ids: List[str] = Field(
        ...,
        description="List of unique business IDs to resolve names for."
    )