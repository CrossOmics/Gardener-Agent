import json
import re
from typing import Dict, Any, List, Optional
from loguru import logger

from langchain.agents import create_agent
from langchain_core.messages import HumanMessage, AIMessage, ToolMessage

from base.base_model_manager import base_llm
from tool.tool_definitions import TOOL_REGISTRY
from base.history_log_call import log_agent_thought
from utils.message_utils import serialize_messages


def _extract_snapshot_id(content: str) -> Optional[str]:
    """Helper to extract snapshot_id from tool output JSON or text."""
    # 1. Try JSON
    try:
        data = json.loads(content)
        if isinstance(data, dict):
            if "snapshot_id" in data:
                return data["snapshot_id"]
            if "final_snapshot_id" in data:
                return data["final_snapshot_id"]
    except (json.JSONDecodeError, TypeError):
        pass

    # 2. Try Regex for text output
    # Matches "Snapshot ID: s_123" or "Final Snapshot ID: s_123"
    # We assume IDs are alphanumeric with underscores
    match = re.search(r"(?:Final\s+)?Snapshot\s+ID:\s*([a-zA-Z0-9_]+)", content, re.IGNORECASE)
    if match:
        return match.group(1)

    return None


async def executor_node(state: Dict[str, Any]) -> Dict[str, Any]:
    """
    Level 2: Executor Node (LLM-Driven).
    Executes the plan step-by-step, maintaining context and recording inputs/outputs.
    """
    session_id = state.get("session_id", "unknown_session")
    project_id = state.get("project_id")
    dataset_id = state.get("dataset_id")
    # Start with the snapshot ID from state, but update it if tools produce new ones
    current_snapshot_id = state.get("snapshot_id")
    filename = state.get("filename")
    plan = state.get("current_plan", {})
    steps = plan.get("phase_2_execution_plan", [])

    if not steps:
        logger.warning("[EXECUTOR] No steps found in plan.")
        return {"next_action": "respond"}

    # Initialize accumulation lists (copy to avoid mutating state directly until return)
    collected_inputs = list(state.get("input") or [])
    collected_outputs = list(state.get("output") or [])

    # Executor memory for this run (stores summaries of previous steps)
    executor_memory: List[str] = []

    # Token counters
    total_input_tokens = 0
    total_output_tokens = 0
    total_tool_calls = 0
    successful_calls = 0

    all_execution_results = []

    logger.info(f"[EXECUTOR] Starting step-by-step execution of {len(steps)} steps...")
    await log_agent_thought(session_id, f"[EXECUTOR] Starting autonomous execution of {len(steps)} steps...")

    for i, step in enumerate(steps):
        step_number = step.get("step", i + 1)
        tool_name = step.get("tool_name")
        reason = step.get("reason")
        suggestion = step.get("suggestion", {})

        logger.info(f"[EXECUTOR] Preparing Step {step_number}: {tool_name}")

        # 1. Validate Tool
        tool = TOOL_REGISTRY.get(tool_name)
        if not tool:
            logger.error(f"[EXECUTOR] Tool '{tool_name}' not found.")
            executor_memory.append(f"Step {step_number} Failed: Tool '{tool_name}' not found.")
            continue

        # 2. Construct Prompt
        memory_context = "\n".join(executor_memory) if executor_memory else "None"

        system_prompt = f"""
You are the **Execution Agent**. You are executing Step {step_number} of a plan.

# GLOBAL CONTEXT
- Project ID: {project_id}
- Dataset ID: {dataset_id}
- Current Snapshot ID: {current_snapshot_id}
- Filename: {filename}

# PREVIOUS STEPS CONTEXT
{memory_context}

# CURRENT STEP INSTRUCTION
- **Tool**: {tool_name}
- **Goal**: {reason}
- **Hints/Context**: {json.dumps(suggestion)}

# RULES
1. **Execute the Tool**: You MUST call the tool '{tool_name}' with appropriate arguments based on the Hints and Context.
2. **Snapshot ID**: If a previous step produced a new Snapshot ID (visible in Context), use it. Otherwise use Current Snapshot ID.
3. **Output**: After the tool runs, provide a brief summary of the result.
"""
        logger.debug(f"[Executor] Step {step_number} Full Prompt:\n{system_prompt}")
        try:
            # 3. Create Agent for this step
            # We bind only the specific tool to reduce confusion and ensure strict adherence to the plan
            agent_app = create_agent(
                model=base_llm,
                tools=[tool],
                system_prompt=system_prompt,
                debug=False
            )

            # 4. Execute
            # We pass a trigger message to start the agent
            result_state = await agent_app.ainvoke({
                "messages": [HumanMessage(content=f"Execute Step {step_number} now.")]
            })

            raw_messages = result_state.get("messages", [])

            # 5. Process Results
            step_input_tokens = 0
            step_output_tokens = 0
            step_tool_calls = 0
            step_success = 0

            step_tool_inputs = []
            step_tool_output_content = ""

            for msg in raw_messages:
                # Metrics
                if isinstance(msg, AIMessage) and getattr(msg, "usage_metadata", None):
                    usage = msg.usage_metadata
                    step_input_tokens += usage.get("input_tokens", 0)
                    step_output_tokens += usage.get("output_tokens", 0)

                # Capture Inputs
                if isinstance(msg, AIMessage) and getattr(msg, "tool_calls", None):
                    step_tool_calls += len(msg.tool_calls)
                    for tc in msg.tool_calls:
                        step_tool_inputs.append(tc.get("args", {}))

                # Capture Outputs & Success
                if isinstance(msg, ToolMessage):
                    content_str = str(msg.content)
                    step_tool_output_content = content_str  # Keep the last tool output

                    # Check success (simplified)
                    if "error" not in content_str.lower() and "fail" not in content_str.lower():
                        step_success += 1

                    # Extract Snapshot ID
                    new_snap_id = _extract_snapshot_id(content_str)
                    if new_snap_id:
                        if new_snap_id not in collected_outputs:
                            collected_outputs.append(new_snap_id)
                        current_snapshot_id = new_snap_id  # Update for next steps
                        logger.info(f"[EXECUTOR] Found new Snapshot ID: {new_snap_id}")

            # Update Accumulators
            total_input_tokens += step_input_tokens
            total_output_tokens += step_output_tokens
            total_tool_calls += step_tool_calls
            successful_calls += step_success

            collected_inputs.extend(step_tool_inputs)

            # Update Memory
            # Get the final agent response or the tool output
            final_msg = raw_messages[-1]
            summary = final_msg.content if isinstance(final_msg.content, str) else "Completed."
            executor_memory.append(f"Step {step_number} ({tool_name}): {summary}")

            # Log
            all_execution_results.extend(raw_messages)
            await log_agent_thought(session_id, f"[EXECUTOR] Step {step_number} completed. Result: {summary}")

        except Exception as e:
            logger.error(f"[EXECUTOR] Step {step_number} failed: {e}")
            executor_memory.append(f"Step {step_number} Failed: {str(e)}")

    # Final Return
    return {
        "execution_results": serialize_messages(all_execution_results),
        "execution_summary": executor_memory,
        "next_action": "reflect",
        "total_input_tokens": state.get("total_input_tokens", 0) + total_input_tokens,
        "total_output_tokens": state.get("total_output_tokens", 0) + total_output_tokens,
        "total_tool_calls": state.get("total_tool_calls", 0) + total_tool_calls,
        "successful_tool_calls": state.get("successful_tool_calls", 0) + successful_calls,
        "input": collected_inputs,
        "output": collected_outputs
    }
