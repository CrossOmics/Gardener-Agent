from typing import Optional, Dict

from pydantic import BaseModel, Field


class ProjectResponse(BaseModel):
    """
    Standard response after project creation.
    Now includes session_id initialized during project setup.
    """
    project_id: str = Field(..., description="The unique business ID of the project")
    project_name: str = Field(..., description="Human-readable project name")
    workspace_path: str = Field(..., description="Local filesystem path for the project")
    session_id: Optional[str] = Field(None, description="The default chat session ID created for the project")
    dataset_id: str = ""
    dataset_name: str = ""

    class Config:
        from_attributes = True


class DatasetEntityMap(BaseModel):
    """
    Nested structure for dataset details and its grouped snapshots.
    """
    dataset_name: str
    # Map structure: { stage_name: { snapshot_id: snapshot_name } }
    snapshots: Dict[str, Dict[str, str]]


class ProjectEntityMapResponse(BaseModel):
    """
    Response DTO containing the full nested entity name map for a project.
    """
    project_id: str
    project_name: str
    datasets: Dict[str, DatasetEntityMap]
