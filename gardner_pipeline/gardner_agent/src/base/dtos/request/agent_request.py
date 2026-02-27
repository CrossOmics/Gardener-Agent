from pydantic import BaseModel
from typing import Optional


class AgentRequest(BaseModel):
    message: str
    snapshot_id: str
    project_id: str
    dataset_id: Optional[str] = None
    session_id: Optional[str] = None
    filename: Optional[str] = None
