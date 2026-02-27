from pydantic import BaseModel
from typing import Optional


class CreateProjectRequest(BaseModel):
    """
    Request body for creating a new empty project workspace.
    """
    project_name: str
    description: Optional[str] = None


class CreateDatasetRequest(BaseModel):
    """
    Request body for importing a raw dataset into an existing project.
    """
    # The ID of the existing project where the dataset will be added
    project_id: str

    # Absolute path to the raw data file on the server's disk
    local_file_path: str

    # Biological metadata (now attached to the dataset, not the project)
    organism: Optional[str] = None
    tissue_type: Optional[str] = None
    description: Optional[str] = None

    # Optional custom name for the dataset (defaults to filename if None)
    dataset_name: Optional[str] = None
