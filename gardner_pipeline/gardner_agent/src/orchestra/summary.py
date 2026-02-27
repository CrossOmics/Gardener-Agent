import re
import json
from typing import Dict, Any, List, Set

from langchain_core.messages import AIMessage
from loguru import logger

from base.base_model_manager import strong_llm
from prompt.role_prompt import SUMMARY_PROMPT
from .state import AgentState
from base.history_log_call import log_agent_final
from base.context.snapshot_context import get_snapshot_info
from base.context.name_context import get_display_names


async def _enrich_and_replace_ids(text_blob: str, id_to_name_map: Dict[str, str]) -> str:
    """Helper to replace all found IDs in a block of text with 'Name (ID)'."""
    if not id_to_name_map:
        return text_blob

    # Sort by length descending to replace longer IDs first (e.g., parent before child)
    sorted_ids = sorted(id_to_name_map.items(), key=lambda item: len(item[0]), reverse=True)

    for entity_id, entity_name in sorted_ids:
        if entity_name and entity_id in text_blob:
            replacement_text = f'{entity_name} ({entity_id})'
            # Use regex with word boundaries to ensure we only replace whole IDs
            id_pattern = re.compile(rf'\b({re.escape(entity_id)})\b')
            text_blob = id_pattern.sub(replacement_text, text_blob)
    return text_blob


async def summary_node(state: AgentState) -> Dict[str, Any]:
    """
    Level 3: Summary Node.
    Generates a final report based on execution results and snapshot details.
    This node enriches the context by replacing system IDs with human-readable names.
    """
    session_id = state.get("session_id", "unknown_session")
    messages = state.get("messages", [])
    user_query = messages[-1].content if messages else "No user query found."

    # Part 1: Consolidate all unstructured text context
    tool_inputs = str(state.get("input", []))
    execution_results = str(state.get("execution_results", []))
    reflection = state.get("reflection_comment", "Analysis completed successfully.")
    generated_snapshot_ids = state.get("output", [])

    # Part 2: Pre-process and enrich structured snapshot details
    snapshot_details_list = []
    all_ids_to_resolve: Set[str] = set()

    # First, gather all IDs from the structured snapshot info
    raw_snapshot_infos = []
    if generated_snapshot_ids:
        logger.info(f"[Summary] Fetching details for {len(generated_snapshot_ids)} snapshots: {generated_snapshot_ids}")
        for snap_id in generated_snapshot_ids:
            try:
                info = await get_snapshot_info(snap_id)
                if isinstance(info, dict):
                    raw_snapshot_infos.append(info)
                    all_ids_to_resolve.add(snap_id)
                    if info.get("parent_snapshot_id"):
                        all_ids_to_resolve.add(info["parent_snapshot_id"])
                    if info.get("dataset_id"):
                        all_ids_to_resolve.add(info["dataset_id"])
                else:
                    snapshot_details_list.append(
                        f"--- Snapshot ID: {snap_id} ---\n(Details unavailable or in unexpected format)")
            except Exception as e:
                logger.warning(f"[Summary] Failed to fetch info for snapshot {snap_id}: {e}")
                snapshot_details_list.append(f"--- Snapshot ID: {snap_id} ---\n(Details unavailable: {str(e)})")

    # Part 3: Gather IDs from all other text sources
    unstructured_text = "\n\n".join([user_query, tool_inputs, execution_results, reflection])
    all_ids_to_resolve.update(re.findall(r'\b((?:s|ds)_[a-zA-Z0-9_.-]+)\b', unstructured_text))

    # Part 4: Bulk resolve all unique IDs found
    id_to_name_map = {}
    if all_ids_to_resolve:
        snapshot_ids = {sid for sid in all_ids_to_resolve if sid.startswith('s_')}
        dataset_ids = {did for did in all_ids_to_resolve if did.startswith('ds_')}
        try:
            if snapshot_ids:
                id_to_name_map.update(await get_display_names(list(snapshot_ids), 'snapshot'))
            if dataset_ids:
                id_to_name_map.update(await get_display_names(list(dataset_ids), 'dataset'))
            logger.info(f"[Summary] Resolved {len(id_to_name_map)} names for {len(all_ids_to_resolve)} unique IDs.")
        except Exception as e:
            logger.warning(f"[Summary] Bulk name resolution failed: {e}")

    # Part 5: Build clean, human-readable snapshot details
    for info in raw_snapshot_infos:
        snap_id = info.get("snapshot_id", "Unknown")
        parent_id = info.get("parent_snapshot_id")
        dataset_id = info.get("dataset_id")

        # Create a cleaner dictionary for the LLM
        display_info = {
            "snapshot": f"{id_to_name_map.get(snap_id, 'Unnamed')} ({snap_id})",
            "parent_snapshot": f"{id_to_name_map.get(parent_id, 'N/A')} ({parent_id})" if parent_id else "None",
            "dataset": f"{id_to_name_map.get(dataset_id, 'N/A')} ({dataset_id})" if dataset_id else "Unknown",
            "created": info.get("create_time"),
            "notes": info.get("user_notes"),
            "parameters": info.get("params_json"),
            "key_outputs": info.get("thumbnail_json", {}).get("summary")
        }
        # Convert to a clean string and append
        snapshot_details_list.append(
            f"--- Snapshot: {id_to_name_map.get(snap_id, snap_id)} ---\n"
            f"{json.dumps(display_info, indent=2)}"
        )

    snapshot_details_text = "\n\n".join(snapshot_details_list)
    if not snapshot_details_text:
        snapshot_details_text = "No new snapshots were generated in this run."

    # Part 6: Enrich all unstructured text as a final sweep
    enriched_user_query = await _enrich_and_replace_ids(user_query, id_to_name_map)
    enriched_tool_inputs = await _enrich_and_replace_ids(tool_inputs, id_to_name_map)
    enriched_execution_results = await _enrich_and_replace_ids(execution_results, id_to_name_map)
    enriched_reflection = await _enrich_and_replace_ids(reflection, id_to_name_map)

    # Part 7: Construct the Full Prompt with enriched context
    full_prompt = (
        f"{SUMMARY_PROMPT}\n\n"
        f"# EXECUTION CONTEXT\n"
        f"## 1. User Request\n{enriched_user_query}\n\n"
        f"## 2. Tool Inputs (Parameters Used)\n{enriched_tool_inputs}\n\n"
        f"## 3. Generated Snapshots (Biological Results)\n{snapshot_details_text}\n\n"
        f"## 4. Execution Logs\n{enriched_execution_results}\n\n"
        f"## 5. Reflection & Quality Assurance\n{enriched_reflection}\n"
    )

    logger.info("[Summary] Generating final biological report with Strong LLM...")
    logger.debug(f" [Summary] Full Prompt {full_prompt}")

    # Part 8. Generate Report
    response = await strong_llm.ainvoke(full_prompt)
    report_content = response.content

    # [DB LOG] 9. Final Response
    await log_agent_final(session_id, text=report_content)
    logger.success("[AGENT_FINAL] Report generated successfully.")

    final_message = AIMessage(content=report_content)
    usage = response.usage_metadata or {}

    # Get the latest snapshot ID from the state to pass it through
    final_snapshot_id = state.get("snapshot_id")

    return {
        "messages": [final_message],
        "next_action": "end",
        "snapshot_id": final_snapshot_id,  # <-- Pass the snapshot_id through
        "total_input_tokens": state.get("total_input_tokens", 0) + usage.get("input_tokens", 0),
        "total_output_tokens": state.get("total_output_tokens", 0) + usage.get("output_tokens", 0)
    }
