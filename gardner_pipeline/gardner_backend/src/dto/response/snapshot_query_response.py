from typing import Dict, Any, Optional
from datetime import datetime
from pydantic import BaseModel


class SnapshotQueryResponse(BaseModel):
    """
    Response DTO for querying a single snapshot's details.
    Separates experiment parameters (input) from biological metrics (output).
    """
    snapshot_id: str
    # Parent snapshot ID, specifying the snapshot on which the current one is based
    parent_snapshot_id: Optional[str] = None
    branch_name: str
    create_time: datetime
    user_notes: Optional[str] = None
    snapshot_name: Optional[str] = None

    # input: The parameters used to run this experiment
    input: Dict[str, Any]

    # output: Key biological metrics (e.g., cell counts, retention rates)
    # File paths (plots/csvs) are filtered out.
    output: Dict[str, Any]
