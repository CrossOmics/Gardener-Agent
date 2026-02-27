from typing import List, Dict, Any
from base.base_call import call_backend


async def get_recent_history(session_id: str, limit: int = 3) -> List[Dict[str, Any]]:
    """
    Retrieves the recent conversation history for a given session.
    Filters specifically for 'user_input' and 'agent_final' message types
    to construct a clean context for the planner.

    Args:
        session_id (str): The unique identifier for the chat session.
        limit (int): The number of conversation turns (user-agent pairs) to retrieve.
                     Defaults to 3.

    Returns:
        List[Dict[str, Any]]: A list of message dictionaries with 'role' ('user' or 'assistant')
                              and 'content'. Returns empty list on error.
    """
    # We request a larger page size (e.g., 20) to ensure we capture enough
    # user/agent messages. The backend filters by type, so this returns
    # 20 messages of type 'user_input' or 'agent_final'.
    endpoint = f"/agent/history/{session_id}/interaction"

    try:
        response = await call_backend("GET", endpoint, payload={"page": 1, "size": 20})
    except Exception as e:
        print(f"[HistoryContext] Failed to fetch history: {e}")
        return []

    if not response or "messages" not in response:
        return []

    raw_messages = response["messages"]
    history = []

    for i in range(len(raw_messages) - 1, -1, -1):
        msg = raw_messages[i]
        if msg is None:
            continue
        ordinal = len(raw_messages) - i
        msg_type = msg.get("type")
        content_payload = msg.get("content", {})

        role = None
        text_content = ""

        if msg_type == "user_input":
            role = "user"
            # Extract text from content.content
            # Structure: {"content": "text"}
            text_content = content_payload.get("content", "")

        elif msg_type == "agent_final":
            role = "assistant"
            # Extract text from content.content list
            # Structure: {"content": [{"type": "text", "text": "..."}]}
            inner_content = content_payload.get("content", [])
            if isinstance(inner_content, list):
                parts = []
                for block in inner_content:
                    if isinstance(block, dict) and block.get("type") == "text":
                        parts.append(block.get("text", ""))
                text_content = "\n".join(parts)
            elif isinstance(inner_content, str):
                text_content = inner_content

        if role and text_content:
            history.append({"ordinal": ordinal, "role": role, "content": text_content})

    # Limit to the requested number of turns.
    # Assuming 1 turn = 1 user message + 1 assistant message (2 messages).
    max_messages = limit * 2

    if max_messages <= 0:
        return []

    if len(history) > max_messages:
        history = history[-max_messages:]

    return history
