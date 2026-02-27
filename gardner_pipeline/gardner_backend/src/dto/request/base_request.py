from typing import Optional
from pydantic import BaseModel, Field


class BaseAnalysisRequest(BaseModel):
    """
    Base Data Transfer Object (DTO) for all analysis requests.
    Contains common identifiers required for context resolution.

    """
    project_id: str = Field(default=None, description="The business ID of the project")
    dataset_id: str = Field(default=None, description="The business ID of the dataset")
    snapshot_id: Optional[str] = Field(default=None, description="The specific snapshot ID to operate on")
