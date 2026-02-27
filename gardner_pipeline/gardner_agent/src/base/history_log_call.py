import copy
from loguru import logger
from typing import Dict, Any, Optional

from base.base_call import call_backend


def _sanitize_data(data: Any) -> Any:
    """
    Iteratively removes 'extras' and 'signature' fields using a stack.
    This avoids recursion limits and is generally efficient.
    """
    try:
        # Create a deep copy to ensure we don't modify the original data used elsewhere
        clean_data = copy.deepcopy(data)
    except Exception as e:
        logger.warning(f"Failed to deepcopy data for sanitization: {e}")
        return data

    # Use a stack for iterative traversal
    stack = [clean_data]

    while stack:
        curr = stack.pop()

        if isinstance(curr, dict):
            # Remove unwanted keys efficiently
            curr.pop("extras", None)
            curr.pop("signature", None)

            # Push nested containers to stack
            for value in curr.values():
                if isinstance(value, (dict, list)):
                    stack.append(value)

        elif isinstance(curr, list):
            # Push items to stack
            for item in curr:
                if isinstance(item, (dict, list)):
                    stack.append(item)

    return clean_data


# Agent History Logging Functions
async def _log_message(session_id: str, message_type: str, data: Dict[str, Any]) -> Optional[str]:
    """
    Internal helper to push log events to the history controller.
    Endpoint: POST /agent/history/{session_id}/message/create
    """
    endpoint = f"/agent/history/{session_id}/message/create"

    # Clean data to remove unnecessary fields like 'extras'
    cleaned_data = _sanitize_data(data)

    # Construct payload matching AgentMessageRequest DTO
    payload = {
        "message_type": message_type,
        "data": cleaned_data
    }

    try:
        # Re-use call_backend for consistency, but handle logging errors gracefully
        # so logging failures don't crash the main agent logic.
        resp = await call_backend("POST", endpoint, payload)
        return resp.get("msg_id")
    except Exception as e:
        logger.warning(f"Failed to log agent event ({message_type}): {e}")
        return None


async def log_user_message(session_id: str, text: str):
    """Logs incoming user query."""
    return await _log_message(session_id, "user_input", {"text": text})


async def log_agent_thought(session_id: str, thought: str):
    """Logs internal reasoning (Chain of Thought)."""
    return await _log_message(session_id, "agent_thought", {"thought": thought})


async def log_tool_call(session_id: str, tool_name: str, arguments: Dict[str, Any]):
    """Logs the intent to call a tool."""
    return await _log_message(session_id, "tool_call", {
        "tool_name": tool_name,
        "arguments": arguments
    })


async def log_tool_result(session_id: str, tool_name: str, result_data: Any, is_error: bool = False):
    """Logs the output observation from a tool."""
    return await _log_message(session_id, "tool_result", {
        "tool_name": tool_name,
        "result_data": result_data,
        "is_error": is_error
    })


async def log_agent_final(session_id: str, text: str, attachments: list = None):
    """Logs the final response displayed to the user."""
    return await _log_message(session_id, "agent_final", {
        "text": text,
        "attachments": attachments or []
    })
