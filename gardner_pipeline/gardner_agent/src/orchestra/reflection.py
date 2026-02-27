import json
from typing import Dict, Any, List

from json_repair import json_repair
from loguru import logger
from langchain_core.messages import HumanMessage, SystemMessage

from base.base_model_manager import base_llm
from base.history_log_call import log_agent_thought
from prompt.role_prompt import REFLECTION_PROMPT


async def reflection_node(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Level 3: Reflection Node.
    Evaluates execution results.
    - If Success: Passes to Summary.
    - If Failure: Generates feedback and passes back to Planner.
    """
    session_id = state.get("session_id", "unknown")
    user_request = state.get("user_request", "")
    execution_results = state.get("execution_results", [])
    execution_summary = state.get("execution_summary", [])
    snapshot_id = state.get("snapshot_id")  # Preserve snapshot_id

    # Extract inputs and outputs from state
    tool_inputs = state.get("input", [])

    generated_snapshots = state.get("output", [])

    retry_count = state.get("retry_count", 0)
    max_retries = 3
    # debug print
    logger.debug(f"tool_inputs: {tool_inputs}, generated_snapshots: {generated_snapshots}")
    logger.info("[REFLECTION] Analyzing execution results...")

    # 1. Check Retry Limit
    if retry_count >= max_retries:
        logger.warning(f"[REFLECTION] Max retries ({max_retries}) reached. Forcing response.")
        return {
            "next_action": "respond",
            "reflection_feedback": f"Max retries reached. The system failed to execute the request after {max_retries} attempts. Please check the logs.",
            "snapshot_id": snapshot_id
        }

    # 2. Prepare Context for LLM
    # Convert execution results (messages) to a readable string log
    logs_text = ""
    if isinstance(execution_results, list):
        for msg in execution_results:
            # Simple serialization for the prompt
            role = getattr(msg, "type", "unknown")
            content = getattr(msg, "content", "")
            logs_text += f"[{role.upper()}]: {content[:500]}...\n"  # Truncate long logs
    else:
        logs_text = str(execution_results)

    # Format the summary list into a readable string
    summary_text = "\n".join(execution_summary) if execution_summary else "No summary available."

    prompt = REFLECTION_PROMPT.format(
        user_request=user_request,
        execution_summary=summary_text,
        execution_logs=logs_text,
        tool_inputs=json.dumps(tool_inputs, indent=2),
        generated_snapshots=json.dumps(generated_snapshots, indent=2)
    )

    logger.debug(f"[Reflection] Full Prompt:\n{prompt}")
    # 3. Call Strong LLM
    try:
        # Use HumanMessage instead of SystemMessage for better compatibility with Gemini
        response = await base_llm.ainvoke([
            HumanMessage(content=prompt)
        ])

        # Handle case where response.content is a list (e.g. Gemini/Anthropic)
        response_content = response.content
        if isinstance(response_content, list):
            full_text = ""
            for block in response_content:
                if isinstance(block, dict) and "text" in block:
                    full_text += block["text"]
                elif isinstance(block, str):
                    full_text += block
            if full_text:
                response_content = full_text

        result = json_repair.loads(response_content)

        is_successful = result.get("is_successful", False)
        next_action = result.get("next_action", "respond")
        advice = result.get("advice", "")
        reasoning = result.get("reasoning", "")

        logger.info(f"[REFLECTION] Decision: {next_action.upper()} | Reason: {reasoning}")
        await log_agent_thought(session_id, f"[REFLECTION] Analysis: {reasoning}. Decision: {next_action}")

        # 4. Handle Actions
        if next_action == "summary":
            # Pass through. The state already has 'input' and 'output' from executor.
            return {
                "next_action": "summary",
                "reflection_feedback": "Execution successful.",
                "snapshot_id": snapshot_id
            }

        elif next_action == "planner":
            # Increment retry count and provide feedback
            return {
                "next_action": "planner",
                "retry_count": retry_count + 1,
                "reflection_feedback": f"PREVIOUS PLAN FAILED. \nREASON: {reasoning}\nADVICE: {advice}\n\nReview the available snapshot info and tool definitions carefully before replanning.",
                "snapshot_id": snapshot_id
            }

        else:
            # Default to respond (Dead end)
            return {
                "next_action": "respond",
                "reflection_feedback": advice,  # This will be used by the responder to explain failure
                "snapshot_id": snapshot_id
            }

    except Exception as e:
        logger.error(f"[REFLECTION] Error during reflection: {e}")
        # Fallback to respond if reflection crashes
        return {
            "next_action": "respond",
            "reflection_feedback": f"System error during reflection: {str(e)}",
            "snapshot_id": snapshot_id
        }
