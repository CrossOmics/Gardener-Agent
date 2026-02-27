from typing import TypedDict, List, Dict, Any, Optional
from langchain_core.messages import BaseMessage



class AgentState(TypedDict):
    messages: List[BaseMessage]
    session_id: Optional[str]
    project_id: Optional[str]
    dataset_id: Optional[str]
    snapshot_id: Optional[str]
    filename: Optional[str]

    current_plan: Optional[Dict[str, Any]]
    next_action: Optional[str]  # 'execute', 'ask_user', 'summary', 'replan'

    execution_results: List[Dict[str, Any]]
    reflection_comment: Optional[str]  # Error message from Reflection

    iteration_count: int

    total_input_tokens: int
    total_output_tokens: int
    total_tool_calls: int
    successful_tool_calls: int

    # input context parameters
    input: List[Any]

    # output results
    output: List[Any]

