from typing import Optional
from pydantic import BaseModel


# Session DTOs
class CreateSessionRequest(BaseModel):
    project_id: str = "p_default"
    initial_name: Optional[str] = "New Bio Analysis"
    agent_name: Optional[str] = "default"
