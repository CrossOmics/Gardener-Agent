from typing import Dict, Any
from langchain_core.messages import AIMessage
from loguru import logger

from base.base_model_manager import base_llm
from .state import AgentState
from prompt.role_prompt import ASK_USER_PROMPT
from base.history_log_call import log_agent_final


async def ask_user_node(state: AgentState) -> Dict[str, Any]:
    """
    Ask User Node - Requests missing information.
    """
    ask_user_llm = base_llm
    session_id = state.get("session_id", "unknown_session")
    user_input = state["messages"][-1].content
    plan = state.get("current_plan", {})

    # Extract context
    context_check = plan.get("phase_1_context_check", {})
    intent_analysis = plan.get("phase_2_intent_analysis", {})

    missing_context = context_check.get("missing_context", [])
    clarification_needed = intent_analysis.get("clarification_needed", "More information needed")
    failure_reason = state.get("reflection_comment", "No error")

    # Format prompt
    prompt = ASK_USER_PROMPT.format(
        missing_context=", ".join(missing_context) if missing_context else "None",
        user_input=user_input,
        clarification_needed=clarification_needed,
        failure_reason=failure_reason
    )

    try:
        response = await ask_user_llm.ainvoke(prompt)
        ask_message = response.content
    except Exception as e:
        logger.error(f"[ASK_USER_ERROR] {e}")
        ask_message = f"I need more information to proceed. Could you provide: {', '.join(missing_context)}?"

    logger.warning(f"[ASK_USER] {ask_message}")

    # Log final response
    await log_agent_final(session_id, text=ask_message)

    final_message = AIMessage(content=ask_message)

    return {
        "messages": [final_message],
        "next_action": "end",
        "iteration_count": 0  # Reset loop counter
    }
