from pydantic import BaseModel
from typing import Dict, Any, Union, List
from pydantic import BaseModel, Field


class AnnotationResultDTO(BaseModel):
    snapshot_id: str
    snapshot_path: str

    # supports a single ID (str) OR a list of clustering ID -> Label (dict)
    cluster_id: Union[str, List, None]
    enrichment_results: Dict[str, Any]
    msg: str


class UpdateAnnotationLabelResponse(BaseModel):
    """
    Standard response after updating cluster annotation labels.
    """
    msg: str = Field(..., description="Status message indicating success or partial success.")
    snapshot_id: str = Field(..., description="The ID of the snapshot that was updated.")
    updated_count: int = Field(..., description="Number of cluster labels that were actually modified.")
    target_file: str = Field(..., description="The specific annotation key (file_name) that was targeted.")
