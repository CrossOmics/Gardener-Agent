from pydantic import BaseModel
from datetime import datetime
from typing import List

class SnapshotAncestorNode(BaseModel):
    """
    DTO for a single node in the snapshot ancestor lineage.
    """
    snapshot_id: str
    snapshot_name: str
    stage_name: str
    create_time: datetime