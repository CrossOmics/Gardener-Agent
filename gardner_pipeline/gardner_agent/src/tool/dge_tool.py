from typing import Literal

from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description=(
        "Run Differential Gene Expression (DGE) analysis to identify marker genes for each cluster. "
        "Use this after clustering to find genes that distinguish different cell groups."
    ),
    parse_docstring=True
)
async def run_dge_analysis(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        groupby: str,
        method: Literal['wilcoxon', 't-test', 'logreg'] = "wilcoxon",
        n_top_genes: int = 25,
        use_raw: bool = True
) -> str:
    """
    Execute Differential Gene Expression (DGE) analysis to identify cluster-specific marker genes.

    # CORE RESPONSIBILITY
    Use this tool to find genes that are significantly upregulated in each cluster compared to all others.
    This is essential for:
    - Understanding what makes each cluster unique
    - Preparing for cell type annotation
    - Identifying potential biomarkers

    # WHEN TO USE
    - After clustering is complete
    - When user asks: "What are the marker genes for cluster X?"
    - When user wants to annotate cell types (DGE is a prerequisite)
    - Before running annotation tools (GSEApy/CellTypist)

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: Ensure clustering has been completed. The `snapshot_id` should point to a clustering result.
    2. **Groupby Parameter**: Usually "leiden" or "louvain" or "cplearn" (the clustering column name). Check snapshot details if unsure.
    3. **Method Selection**:
       - `wilcoxon`: Default, robust for most cases
       - `t-test`: Faster but assumes normal distribution
       - `logreg`: Good for binary comparisons
    4. **Top Genes**: Default 25 is usually sufficient. Increase if user wants more detailed analysis.

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a clustering snapshot (e.g., 's_cluster_leiden_...').
    - **`groupby` must exist**: The column name must exist in the AnnData object. (Options: ('leiden', 'louvain', 'cplearn').

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from clustering result (e.g., 's_cluster_leiden_...'). REQUIRED.
        groupby: Column name in adata.obs to group by (Options: ('leiden', 'louvain', 'cplearn').
        method: Statistical test method. Options: "wilcoxon" (default, robust), "t-test" (faster), "logreg" (binary).
        n_top_genes: Number of top marker genes to identify per cluster (5-100). Default: 25.
        use_raw: Whether to use adata.raw for calculation (recommended: True). Default: True.

    Returns:
        Summary of DGE analysis results including:
        - Number of clusters analyzed
        - Top marker genes per cluster (preview)
        - Paths to visualization plots (rank plot, dot plot, heatmap, violin plot)
        - CSV file path with full statistical results
        - New snapshot ID
    """
    try:
        data = await call_backend("POST", "/dge/new", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "groupby": groupby,
            "method": method,
            "n_top_genes": n_top_genes,
            "use_raw": use_raw
        })

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        top_markers = data.get("top_markers", {})
        csv_path = data.get("csv_path", "")
        
        # Format top markers preview
        markers_summary = []
        for cluster_id, genes in list(top_markers.items())[:5]:  # Show first 5 clusters
            genes_str = ", ".join(genes[:5])  # Show first 5 genes per cluster
            markers_summary.append(f"  Cluster {cluster_id}: {genes_str}")

        return (
            f"DGE Analysis Complete!\n"
            f"\n"
            f"- Method: {method}\n"
            f"- Grouped by: {groupby}\n"
            f"- Clusters analyzed: {len(top_markers)}\n"
            f"- Top genes per cluster: {n_top_genes}\n"
            f"\n"
            f"Top Marker Genes Preview:\n"
            f"{chr(10).join(markers_summary)}\n"
            f"\n"
            f"Results:\n"
            f"  - Full CSV: {csv_path}\n"
            f"  - Rank plot: {data.get('rank_genes_plot', 'N/A')}\n"
            f"  - Dot plot: {data.get('dotplot', 'N/A')}\n"
            f"  - Violin plot: {data.get('violin', 'N/A')}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review marker genes to understand cluster identity\n"
            f"  2. Run annotation tools (GSEApy or CellTypist) to assign cell types\n"
            f"  3. Use the CSV file for detailed statistical analysis"
        )
    except Exception as e:
        return f"DGE analysis failed: {str(e)}"