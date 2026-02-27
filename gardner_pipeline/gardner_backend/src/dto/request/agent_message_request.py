from typing import Dict, Any

from pydantic import BaseModel, Field


# Agent Message DTO
class AgentMessageRequest(BaseModel):
    """
    Unified payload for logging any agent event.
    The 'payload' structure depends on 'message_type'.
    """
    message_type: str = Field(
        ...,
        description="Event type: user_input, agent_thought, tool_call, tool_result, agent_final"
    )
    data: Dict[str, Any] = Field(
        ...,
        description="Dynamic payload matching the message type requirements."
    )
