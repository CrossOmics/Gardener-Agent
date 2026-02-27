from pydantic import BaseModel
from typing import Optional


class AgentResponse(BaseModel):
    reply: str
    final_snapshot_id: Optional[str] = None
