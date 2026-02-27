from typing import List, Union

from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description="Execute clustering (Leiden). Use this to merge fragmented clusters (lower resolution) "
                "or split mixed clusters (higher resolution).",
    parse_docstring=True
)
async def run_clustering(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        method: str = "leiden",
        resolution: float = 0.5
) -> str:
    """
    Execute clustering analysis with a specific resolution to reorganize cell groups.

    # CORE RESPONSIBILITY
    Use this tool when the researcher is unsatisfied with the current cell grouping.
    - If they say "clusters are too fragmented" or "merge groups" -> **Decrease** resolution.
    - If they say "clusters are mixed" or "split groups" -> **Increase** resolution.

    # STRATEGY & BEST PRACTICES
    1. **Check Context First**: Before calling this, you MUST know the *current* resolution.
       Call `get_snapshot_details` first if you don't know it.
    2. **Tuning Logic**:
       - To Merge: New Resolution = Current Resolution - 0.2 (approx).
       - To Split: New Resolution = Current Resolution + 0.2 (approx).

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: You must strictly use the `project_id` from the conversation context.
      Do not invent or guess this ID. If missing, ask the user.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset to be clustered.
        snapshot_id: Snapshot ID to cluster from (typically the 'Neighborhood Graph' snapshot).
        method: Clustering algorithm. Options['leiden', 'louvain', 'cplearn']. Follow the user intented method.
        resolution: Granularity parameter. Higher (e.g., 1.0) = more clusters. Lower (e.g., 0.3) = fewer clusters.
    """
    try:
        # Implementation remains the same, calling the backend service via HTTP
        data = await call_backend("POST", "/clustering/create", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "method": method,
            "resolution": resolution
        })

        # Extract summary for the agent's observation
        n_clusters = list(data.get("clusters_summary", {}).values())[0]
        return (
            f"Clustering complete.\n"
            f"- Method: {method}\n"
            f"- New Resolution: {resolution}\n"
            f"- Clusters Found: {n_clusters}\n"
            f"- New Snapshot ID: {data.get('snapshot_id')}\n"
            f"Observation: The data has been re-clustered. ask the user to verify if the separation is better."
        )
    except Exception as e:
        return f"Clustering failed: {str(e)}"


@tool(
    description="Merge specific clusters manually into a single group. Use when the researcher explicitly identifies clusters to combine (e.g., 'combine cluster 0 and 1').",
    parse_docstring=True
)
async def merge_clusters(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        method: str,
        clusters_to_merge: Union[List[str], List[int], List[float]],
        new_label: str = "Merged_Group",
) -> str:
    """
    Manually merge multiple specific clusters into a single new cluster label.

    # CORE RESPONSIBILITY
    Use this tool for precise corrections when the researcher wants to combine specific groups without re-running the entire clustering algorithm.

    # PARAMETER GUIDELINES
    1. **snapshot_id**: This is CRITICAL. You must look back in the conversation history to find the **Snapshot ID of the CLUSTERING result** you are currently viewing.
       - It usually looks like `s_cluster_leiden_...` or `s_cluster_louvain_...`.
       - Do NOT use the dataset_id here. Do NOT invent a new ID. Find the existing one.
    2. **clusters_to_merge**: Extract the exact IDs the user mentioned (e.g., [0, 1] or ["3", "5"]). Both numbers and strings are accepted.
    3. **new_label**: If the user didn't provide a name (e.g., "T-cells"), you can leave this blank, and it will default to "Merged_Group".

    Args:
        project_id: Unique business ID of the project. REQUIRED.
        dataset_id: Unique ID of the dataset.
        method: The clustering key in the data to modify. Options['leiden', 'louvain', 'cplearn']. Follow the user intented method.
        snapshot_id: The Source Snapshot ID of the *existing* clustering result to be modified (e.g. 's_cluster_...'). .
        clusters_to_merge: A list of cluster IDs to combine.
        new_label: The cluster name for the new merged group. Optional parameter

    """
    try:
        # Convert any numbers to strings to ensure backend compatibility
        # (Though backend handles it, doing it here is safer for JSON serialization)
        safe_clusters = [str(c) for c in clusters_to_merge]

        data = await call_backend("POST", "/clustering/merge", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "method": method,
            "clusters_to_merge": safe_clusters,
            "new_label": new_label
        })

        # Extract info for agent observation
        n_clusters = list(data.get("clusters_summary", {}).values())[0]
        new_snap = data.get("snapshot_id")

        return (
            f"Clusters merged successfully.\n"
            f"- Merged IDs: {', '.join(safe_clusters)} -> '{new_label}'\n"
            f"- Total Clusters Remaining: {n_clusters}\n"
            f"- New Snapshot ID: {new_snap}\n"
            f"Observation: The clusters have been merged into a new snapshot. The UMAP has been updated to show '{new_label}'."
        )
    except Exception as e:
        return f"Merge operation failed: {str(e)}"


@tool(
    description="Sub-cluster specific groups. Use this to isolate one or more clusters "
                "and re-run the clustering algorithm on just that subset to split them further.",
    parse_docstring=True
)
async def run_sub_clustering(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        method: str,
        target_clusters: Union[List[str], List[int]],
        resolution: float = 0.5,
        overwrite: bool = False
) -> str:
    """
    Execute sub-clustering on specific target clusters to reveal finer substructures.

    # CORE RESPONSIBILITY
    Use this tool when the researcher identifies that specific clusters are too "coarse" or contain mixed cell types.
    - Example: "Cluster 1 looks like it contains both CD4 and CD8 T-cells. Please sub-cluster it."
    - This creates a hierarchy: Cluster '1' becomes '1-0', '1-1', etc.

    # STRATEGY & BEST PRACTICES
    1. **Target Identification**: You must identify exactly which cluster IDs the user wants to drill down into (e.g., ["1", "3"]).
    2. **Resolution**:
       - Default is 0.5.
       - Higher resolution (e.g., 1.0) forces more splits within the sub-population.
    3. **Snapshot Context**: This operation requires the *result* of a previous clustering run (the `snapshot_id` where the target clusters currently exist).

    Args:
        project_id: Unique business ID of the project. REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: The source snapshot ID containing the clusters to be split. REQUIRED.
        method: Clustering algorithm to apply to the subset. Have value: ['leiden', 'louvain', 'cplearn']. Follow the user intented method.
        target_clusters: A list of cluster IDs to isolate and re-cluster (e.g., ["1", "5"]).
        resolution: Granularity for the sub-clustering. Defaults to 0.5.
        overwrite: If True, updates the existing snapshot in-place. If False, creates a new analysis branch. Defaults to False.
    """
    try:
        # 1. Sanitize inputs: Ensure all cluster IDs are strings for the backend
        safe_targets = [str(c) for c in target_clusters]

        # 2. Call the backend API
        data = await call_backend("POST", "/clustering/sub", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "target_clusters": safe_targets,
            "method": method,
            "resolution": resolution,
            "overwrite": overwrite
        })

        # 3. Format the observation for the agent
        # The backend returns a summary of the *new* total clusters or the specific sub-clusters
        n_clusters = list(data.get("clusters_summary", {}).values())[0]
        new_snap = data.get("snapshot_id")
        action_type = "Overwritten" if overwrite else "Created New Branch"

        return (
            f"Sub-clustering complete ({action_type}).\n"
            f"- Targeted Clusters: {safe_targets}\n"
            f"- Method: {method} (Res: {resolution})\n"
            f"- Total Clusters in Dataset: {n_clusters}\n"
            f"- Snapshot ID: {new_snap}\n"
            f"Observation: The specified clusters have been split. The grouping labels now follow the 'Parent-Child' format (e.g., '1-0')."
        )

    except Exception as e:
        return f"Sub-clustering failed: {str(e)}"
