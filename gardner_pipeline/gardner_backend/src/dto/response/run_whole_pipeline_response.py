from pydantic import BaseModel, Field
from typing import List, Dict


class RunWholePipelineResponse(BaseModel):
    """
    Response payload returning the execution status and generated snapshots.
    """
    pipeline_run_id: str = Field(..., description="Unique ID for this pipeline execution")
    final_snapshot_id: str = Field(..., description="Snapshot ID of the final result")
    status: str = Field("SUCCESS", description="Execution status")
    msg: str = Field("Whole pipeline executed successfully.", description="Result message")

    # List of completed step names
    steps_completed: List[str] = Field(
        ...,
        example=["QC", "HVG", "PCA", "Neighbors", "Clustering", "DEG", "Annotation"]
    )

    # Dictionary mapping step names to their snapshot IDs
    snapshots: Dict[str, str] = Field(
        ...,
        example={"qc": "s_qc_1", "clustering": "s_cluster_2"}
    )
