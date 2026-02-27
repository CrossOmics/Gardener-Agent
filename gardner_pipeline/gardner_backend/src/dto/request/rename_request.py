from typing import Literal

from pydantic import BaseModel, Field


class RenameRequest(BaseModel):
    """
    Data Transfer Object for renaming a project, dataset, or snapshot display name.
    """
    id_type: Literal['project', 'dataset', 'snapshot'] = Field(
        ..., description="The type of entity being renamed."
    )
    current_id: str = Field(
        ..., description="The unique system identifier of the entity."
    )
    new_name: str = Field(
        ..., min_length=1, description="The new user-friendly display name."
    )
