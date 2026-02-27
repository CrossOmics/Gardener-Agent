from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description=(
            "Retrieve metadata and parameters of an existing analysis snapshot. "
            "Use this tool to understand how a snapshot was generated "
            "(e.g. clustering method and resolution) before running further analysis."
    ),
    parse_docstring=True
)
async def get_snapshot_details(snapshot_id: str) -> str:
    """
    Fetch metadata and parameters for a specific analysis snapshot.

    # CORE RESPONSIBILITY
    Use this tool when you need to "look back" at the history of the data.
    It is ESSENTIAL for grounding your next decision in existing parameters.
    - If the user says "what resolution was used?" -> Call this.
    - If you need to re-cluster but don't know the current resolution -> Call this.

    # STRATEGY & BEST PRACTICES
    1. Grounding: Always call this before proposing changes to resolution or methods.
       This prevents redundant work or illogical parameter jumps.
    2. Context Building: Use the returned 'Method' and 'Resolution' to explain
       the current status to the researcher.

    # CONSTRAINTS
    - Requires a valid `snapshot_id`. If the ID is 'Unknown' or missing, you must
      find it via `query_latest_stage_snapshot` first.

    Args:
        snapshot_id: Unique identifier of the analysis snapshot to inspect (e.g., 's_cluster_...'). REQUIRED.
    """
    try:
        data = await call_backend("GET", f"/snapshots/retrieve/{snapshot_id}")

        params = data.get("params_json", {})
        method = params.get("method", "unknown")
        resolution = params.get("resolution", "unknown")
        full_resolutions = params.get("resolutions", [])

        return (
            "Snapshot Details retrieved successfully:\n"
            f"- Project ID: {data.get('project_id')}\n"
            f"- Dataset ID: {data.get('dataset_id')}\n"
            f"- Method: {method}\n"
            f"- Current Resolution: {resolution}\n"
            f"- Evaluated Resolutions: {full_resolutions}\n"
            f"Observation: Now you know the baseline. Proceed with parameter adjustments if requested."
        )
    except Exception as e:
        return f"Failed to retrieve snapshot details for '{snapshot_id}': {str(e)}"


@tool(
    description="Delete a specific analysis snapshot and its associated files from the system.",
    parse_docstring=True
)
async def delete_snapshot(
        snapshot_id: str
) -> str:
    """
    Permanently delete a snapshot record and its physical data files.

    # CORE RESPONSIBILITY
    Use this tool only when the user explicitly asks to remove a specific analysis step
    to save space or clean up the history.

    # CONSTRAINTS
    - This action is IRREVERSIBLE.
    - You must have a specific `snapshot_id` (e.g., 's_cluster_leiden_...').

    Args:
        snapshot_id: The unique business ID of the snapshot to be deleted. REQUIRED.
    """
    try:
        # Communicate with the backend via DELETE method
        data = await call_backend("DELETE", f"/snapshots/delete/{snapshot_id}")

        return f"Snapshot '{snapshot_id}' has been successfully deleted. Result: {data.get('msg')}"
    except Exception as e:
        return f"Failed to delete snapshot: {str(e)}"


@tool(
    description="Retrieve the latest snapshot for a specific dataset and analysis stage (e.g., 'Preprocessing').",
    parse_docstring=True
)
async def query_latest_stage_snapshot(
        dataset_id: str,
        branch_name: str
) -> str:
    """
    Fetch the most recent snapshot ID and its metrics for a specific analysis stage.

    # CORE RESPONSIBILITY
    Use this tool when you need to find the "current state" of an analysis branch without
    having a specific Snapshot ID from the context.

    # STRATEGY
    - Valid branch_names include: 'Preprocessing', 'Clustering', 'DGE', 'Annotation'.
    - Use this to verify if a stage has already been completed.

    Args:
        dataset_id: Unique business ID of the dataset. REQUIRED.
        branch_name: The analysis stage name (e.g., 'Clustering'). REQUIRED.
    """
    try:
        # Querying the latest state using URL parameters
        data = await call_backend("GET", "/snapshots/stage/query/", {
            "dataset_id": dataset_id,
            "branch_name": branch_name
        })

        snapshot_id = data.get("snapshot_id")
        user_notes = data.get("user_notes", "N/A")

        return (
            f"Latest snapshot for stage '{branch_name}' found.\n"
            f"- Snapshot ID: {snapshot_id}\n"
            f"- Created At: {data.get('create_time')}\n"
            f"- User Notes: {user_notes}\n"
            f"Observation: This is the most current data for {branch_name}."
        )
    except Exception as e:
        return f"Failed to query latest snapshot for stage '{branch_name}': {str(e)}"
