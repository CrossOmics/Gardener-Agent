from typing import Dict, Any, List

from langchain_core.messages import HumanMessage, AIMessage, ToolMessage


def serialize_messages(messages: List[Any]) -> List[Dict[str, Any]]:
    """
    Helper to convert LangChain messages to JSON-serializable dicts.
    Extracts text content and tool call details.
    """
    serialized = []
    for m in messages:
        if isinstance(m, HumanMessage):
            serialized.append({
                "role": "user",
                "content": m.content
            })
        elif isinstance(m, AIMessage):
            # Capture tool calls if present
            tool_calls = m.tool_calls if hasattr(m, "tool_calls") else []
            serialized.append({
                "role": "assistant",
                "content": m.content,
                "tool_calls": tool_calls
            })
        elif isinstance(m, ToolMessage):
            serialized.append({
                "role": "tool",
                "name": m.name,
                "content": m.content,
                "tool_call_id": m.tool_call_id
            })
    return serialized
