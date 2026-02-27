
from typing import Dict, Any
from langchain_core.messages import AIMessage
from loguru import logger
from .state import AgentState
from base.history_log_call import log_agent_final


async def respond_node(state: AgentState) -> Dict[str, Any]:
    """
    Direct Response Node - Handles informational queries without tool execution.
    """
    session_id = state.get("session_id", "unknown_session")
    plan = state.get("current_plan", {})
    direct_response = plan.get("direct_response", "I'm here to help with your analysis!")

    logger.info(f"[Quick Respond] Sending direct response: {direct_response[:100]}...")

    # Log to database
    await log_agent_final(session_id, text=direct_response)

    # Return as AI message
    final_message = AIMessage(content=direct_response)
    # Debug Log
    logger.debug(f"[Quick Respond] {final_message}")
    return {
        "messages": [final_message],
        "next_action": "end"
    }
