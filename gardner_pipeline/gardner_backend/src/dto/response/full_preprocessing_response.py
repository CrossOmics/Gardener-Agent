from pydantic import BaseModel
from typing import List


class FullPreprocessingResponse(BaseModel):
    """
    Response DTO for one-shot preprocessing pipeline execution.
    """
    project_id: str
    dataset_id: str

    # Standard Snapshot Identification
    snapshot_id: str
    snapshot_path: str

    # Pipeline execution summary
    executed_steps: List[str]
    skipped_steps: List[str]

    # Semantic state of the final AnnData
    final_stage: str

    message: str = "Preprocessing pipeline completed successfully."
