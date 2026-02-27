from typing import Dict, Any
from json_repair import json_repair
from langchain_core.messages import SystemMessage

from loguru import logger

from base.base_model_manager import base_llm, strong_llm
from base.context.snapshot_context import get_snapshot_info, get_snapshot_ancestors_full, get_snapshot_peers
from base.context.history_message_context import get_recent_history
from .state import AgentState
from prompt.role_prompt import PLANNER_ROLE_PROMPT
from tool.tool_definitions import get_tool_definitions_str
from base.history_log_call import log_user_message, log_agent_thought


async def planner_node(state: AgentState) -> Dict[str, Any]:
    """
    Level 1: Planner Node with THREE-WAY ROUTING
    Routes to: respond / ask_user / execute
    """
    # Use Strong LLM for complex planning
    planner_llm = strong_llm
    session_id = state.get("session_id", "unknown_session")
    current_iter = state.get("iteration_count", 0) + 1

    # Loop Protection
    if current_iter > 3:
        logger.warning(f"[LOOP DETECTED] Max retries exceeded. Routing to ask_user.")
        return {
            "iteration_count": 0,
            "next_action": "ask_user",
            "reflection_comment": f"Failed after {current_iter - 1} attempts. Last error: {state.get('reflection_comment')}"
        }

    # Prepare Context
    user_input = state["messages"][-1].content
    # Fetch Recent Conversation History
    history_context_str = ""
    try:
        # Fetch last 5 turns (user + agent)
        recent_msgs = await get_recent_history(session_id, limit=5)
        if recent_msgs is not None:
            history_context_str = f"""
    ## 5. Recent Conversation History
    (Most recent messages, use this to understand context, follow-ups, or 'wrap up' requests)
    {recent_msgs}
    end of history context. ordinal 1 -> n, earlier to latest
    """
            logger.debug(f"[Planner] Loaded {len(recent_msgs)} recent history messages.")
    except Exception as e:
        logger.warning(f"[Planner] Failed to load history context: {e}")

    if current_iter == 1:
        await log_user_message(session_id, user_input)
        logger.info(f"[USER_INPUT] {user_input}")

    tool_definitions_str = get_tool_definitions_str()

    # Build Failure Context for Replanning
    failure_context = ""
    previous_error = state.get("reflection_comment")
    if previous_error and current_iter > 1:
        failure_context = f"""
## PREVIOUS ATTEMPT FAILED
**Error:** {previous_error}

**Recovery Instructions:**
1. Fix the issue described above
2. DO NOT repeat the same plan
"""

    # Dynamic Snapshot Context Construction
    snapshot_id = state.get('snapshot_id')
    snapshot_context_str = ""

    if snapshot_id and snapshot_id != "Unknown":
        try:
            # Fetch context concurrently or sequentially
            # Using sequential for simplicity, but could be gathered
            current_info = await get_snapshot_info(snapshot_id)
            ancestors = await get_snapshot_ancestors_full(snapshot_id)
            peers = await get_snapshot_peers(snapshot_id)

            snapshot_context_str = f"""
# DYNAMIC SNAPSHOT CONTEXT
## Current Snapshot
{current_info}

## Snapshot Ancestors (Lineage)
{ancestors}

## Snapshot Peers (Same Branch)
{peers}
"""
            logger.info(f"[Planner] Loaded context for snapshot {snapshot_id}")
        except Exception as e:
            logger.warning(f"[Planner] Failed to load snapshot context: {e}")

    # Format System Prompt with Dynamic Context
    # We append the snapshot context to the end of the prompt
    base_prompt = PLANNER_ROLE_PROMPT.format(
        tool_definitions=tool_definitions_str,
        project_id=state.get('project_id', 'Unknown'),
        dataset_id=state.get('dataset_id', 'Unknown'),
        snapshot_id=state.get('snapshot_id', 'Unknown'),
        filename=state.get('filename', 'Unknown'),
        failure_context=failure_context
    )

    # Append both snapshot context and history context
    system_prompt = snapshot_context_str + history_context_str + base_prompt

    full_system_msg = SystemMessage(content=system_prompt)
    messages = [full_system_msg] + state["messages"]

    # Debug: Print the full prompt before invoking LLM
    logger.debug(f"[Planner] Full Prompt:\n{system_prompt}")

    # Invoke LLM
    try:
        response = await planner_llm.ainvoke(messages)
        usage = response.usage_metadata
        input_tokens = usage.get("input_tokens", 0)
        output_tokens = usage.get("output_tokens", 0)
        response_content = response.content

        # Handle case where response.content is a list (e.g. Gemini/Anthropic)
        if isinstance(response_content, list):
            full_text = ""
            for block in response_content:
                if isinstance(block, dict) and "text" in block:
                    full_text += block["text"]
                elif isinstance(block, str):
                    full_text += block
            if full_text:
                response_content = full_text

        logger.debug(response_content)
        plan_json = json_repair.loads(response_content)

        # DEBUG Log
        # logger.debug(f"[Planner] {json.dumps(plan_json, indent=2)}")
        phase_1_intent_analysis = plan_json.get("phase_1_intent_analysis")
        # Extract Key Fields
        action = phase_1_intent_analysis.get("action", "informational")  # Default to execute if unclear
        direct_resp = plan_json.get("direct_response")
        user_goal = plan_json.get('phase_1_intent_analysis', {}).get('user_goal', 'Unknown goal')

        # Logging
        if direct_resp:
            await log_agent_thought(session_id, f"[RESPOND MODE] Direct answer: {direct_resp[:100]}...")
            logger.info(f"[PLANNER] → RESPOND")
        elif action == "ask_user":
            clarification = plan_json.get('phase_2_intent_analysis', {}).get('clarification_needed', 'Unknown')
            await log_agent_thought(session_id, f"[ASK_USER MODE] Need: {clarification}")
            logger.warning(f"[PLANNER] → ASK_USER: {clarification}")
        else:
            await log_agent_thought(session_id, f"[EXECUTE MODE] Plan (Iter {current_iter}): {user_goal}")
            logger.info(f"[PLANNER] → EXECUTE: {user_goal}")

        # THREE-WAY ROUTING
        return {
            "current_plan": plan_json,
            "next_action": action,  # "respond" / "execute"
            "iteration_count": current_iter,
            "total_input_tokens": state.get("total_input_tokens", 0) + input_tokens,
            "total_output_tokens": state.get("total_output_tokens", 0) + output_tokens
        }

    except Exception as e:
        logger.error(f"[PLANNING_ERROR] {str(e)}")
        return {
            "next_action": "ask_user",
            "reflection_comment": f"Internal planning error: {str(e)}"
        }
