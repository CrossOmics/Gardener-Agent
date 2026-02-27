from typing import List, Optional, Dict
from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description="Run a full annotation pipeline using GSEApy (Enrichr) and/or CellTypist to identify cell types.",
    parse_docstring=True
)
async def run_full_annotation_pipeline(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        target_cluster_col: str,
        categories: List[str] = ["CellMarker_2024"],
        model_names: List[str] = ["Immune_All_Low.pkl"],
        top_n_genes: int = 100,
        majority_voting: bool = True,

) -> str:
    """
    Consolidated annotation pipeline that identifies cell types for cell clusters.

    # CORE RESPONSIBILITY
    Use this tool when the researcher wants to assign biological names (cell types) to the clusters.
    - If they want biological enrichment -> Use GSEApy (pass 'categories').
    - If they want automated model-based classification -> Use CellTypist (pass 'model_names').

    # STRATEGY & BEST PRACTICES
    1. Check Input Snapshot: This tool typically runs on a 'DGE' snapshot (containing marker genes)
       or a 'Clustering' snapshot.
    2. Default Choices:
       - Use 'CellMarker_2024' for GSEApy if the user doesn't specify libraries.
       - Use 'Immune_All_Low.pkl' for CellTypist if working with blood/immune data.
    3. The model_names and categories fields may be inferred from the filename parameter, typically by removing the .png suffix.

    Args:
        project_id: Unique business ID of the project. REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Source Snapshot ID (DGE or Clustering).
        categories: List of Enrichr libraries (options, ['CellMarker_2024', 'NCI-Nature_2016', 'GO_Biological_Process_2023']). Pass empty to skip GSEApy.
        model_names: List of CellTypist models. (options., ['Immune_All_Low.pkl', 'Adult_Human_MTG.pkl', 'Human_Lung_Atlas.pkl']) Pass empty to skip CellTypist.
        top_n_genes: Number of top marker genes to use for GSEApy.
        majority_voting: Whether to refine CellTypist results via cluster-based voting.
        target_cluster_col: The obs column containing cluster IDs. Options: ('leiden', 'louvain', 'cplearn').
    """
    try:
        data = await call_backend("POST", "/annotation/full", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "categories": categories,
            "model_names": model_names,
            "top_n_genes": top_n_genes,
            "majority_voting": majority_voting,
            "target_cluster_col": target_cluster_col
        })

        new_snap = data.get("snapshot_id")
        return (
            f"Annotation pipeline complete.\n"
            f"- GSEApy categories run: {categories}\n"
            f"- CellTypist models run: {model_names}\n"
            f"- New Snapshot ID: {new_snap}\n"
            f"Observation: Cell types have been predicted. You should now display the results "
            f"and ask the user if the labels look correct."
        )
    except Exception as e:
        return f"Annotation pipeline failed: {str(e)}"


@tool(
    description="Manually update/correct predicted cell type labels in an annotation snapshot.",
    parse_docstring=True
)
async def update_annotation_labels(
        snapshot_id: str,
        file_name: str,
        updated_annotation: Dict[str, str]
) -> str:
    """
    Updates the biological labels for clusters based on researcher feedback.

    # CORE RESPONSIBILITY
    Use this tool when the researcher disagrees with the automated predictions and says
    things like: "Cluster 0 is actually T-Cells, not B-Cells".

    # STRATEGY
    1. Identification: You must know the 'file_name' (the key in the metadata) that you are updating.
    2. Map Construction: Provide a dictionary mapping the cluster ID to the new string label.

    Args:
        snapshot_id: The ID of the existing annotation snapshot. REQUIRED.
        file_name: The target annotation key (e.g., 'gseapy_CellMarker_2024' or 'celltypist_majority_voting'). This could
        be inferred from the input filename parameter, typically by removing the .png suffix.

        updated_annotation: A map of cluster IDs to new labels. e.g., {"0": "CD4+ T-cell", "1": "Monocyte"}.
    """
    try:
        data = await call_backend("PATCH", "/annotation/labels/update", {
            "snapshot_id": snapshot_id,
            "file_name": file_name,
            "updated_annotation": updated_annotation
        })

        return (
            f"Labels updated successfully in snapshot '{data.get('snapshot_id')}'.\n"
            f"- Target: {data.get('target_file')}\n"
            f"- Count: {data.get('updated_count')} cluster(s) modified.\n"
            f"Observation: The manual corrections have been applied to the snapshot metadata."
        )
    except Exception as e:
        return f"Label update failed: {str(e)}"


@tool(
    description=(
            "Run GSEApy functional enrichment analysis for ALL clusters in batch. "
            "Use this to automatically annotate all clusters at once using marker genes from DGE analysis."
    ),
    parse_docstring=True
)
async def run_all_gseapy_annotation(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        top_n_genes: int = 100,
        categories: Optional[List[str]] = None
) -> str:
    """
    Execute batch GSEApy annotation for all clusters using marker genes.

    # CORE RESPONSIBILITY
    This tool automatically annotates ALL clusters by:
    1. Loading DGE analysis results (marker genes for each cluster)
    2. Querying Enrichr databases (e.g., CellMarker) for each cluster
    3. Identifying the most likely cell type for each cluster
    4. Creating a summary mapping: Cluster ID -> Cell Type & Confidence

    # WHEN TO USE
    - After DGE analysis is complete
    - When user wants to annotate all clusters at once
    - When user asks: "What cell types are these clusters?"
    - Before detailed single-cluster analysis

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: DGE analysis must be completed first. The `snapshot_id` should point to a DGE result snapshot.
    2. **Categories**: User-friendly database names. Options include:
       - "cellmarker": Cell type markers (most common)
       - "functional": GO and KEGG pathways
       - "msigdb": MSigDB hallmark gene sets
       - If None, uses all available databases
    3. **Top N Genes**: Number of marker genes to use per cluster (default 100). More genes = more comprehensive but slower.

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a DGE analysis snapshot (e.g., 's_dge_wilcoxon_...').
    - **DGE Required**: This tool requires DGE analysis to be completed first.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from DGE analysis result (e.g., 's_dge_wilcoxon_...'). REQUIRED.
        top_n_genes: Number of top marker genes to use per cluster for enrichment (default: 100).
        categories: List of database categories to query. Examples: ["CellMarker_2024", "BioPlanet_2019", ""].

    Returns:
        Summary of batch annotation results including:
        - Number of clusters annotated
        - Cell type predictions for each cluster
        - Confidence scores
        - New snapshot ID
    """
    try:
        payload = {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "top_n_genes": top_n_genes,
            "cluster_id": "ALL"
        }
        if categories:
            payload["categories"] = categories

        data = await call_backend("POST", "/annotation/gseapy/all", payload)

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        enrichment_results = data.get("enrichment_results", {})
        msg = data.get("msg", "")

        # Format cluster annotations summary
        annotations_summary = []
        if isinstance(enrichment_results, dict):
            for cluster_id, result in list(enrichment_results.items())[:10]:  # Show first 10 clusters
                if isinstance(result, dict):
                    cell_type = result.get("cell_type", "Unknown")
                    confidence = result.get("confidence", "N/A")
                    annotations_summary.append(f"  Cluster {cluster_id}: {cell_type} (confidence: {confidence})")

        return (
            f"Batch GSEApy Annotation Complete!\n"
            f"\n"
            f"- Clusters annotated: {len(enrichment_results) if isinstance(enrichment_results, dict) else 'N/A'}\n"
            f"- Top genes per cluster: {top_n_genes}\n"
            f"- Databases used: {', '.join(categories) if categories else 'All available'}\n"
            f"\n"
            f"Cell Type Predictions:\n"
            f"{chr(10).join(annotations_summary) if annotations_summary else '  (Results processing...)'}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Message: {msg}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review cell type predictions for each cluster\n"
            f"  2. Use single-cluster annotation for detailed analysis of specific clusters\n"
            f"  3. Consider running CellTypist for alternative annotation method"
        )
    except Exception as e:
        return f"Batch GSEApy annotation failed: {str(e)}"


@tool(
    description=(
            "Run GSEApy functional enrichment analysis for a SINGLE specific cluster. "
            "Use this to get detailed annotation for one cluster of interest."
    ),
    parse_docstring=True
)
async def run_single_gseapy_annotation(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        cluster_id: str,
        top_n_genes: int = 100,
        categories: Optional[List[str]] = None
) -> str:
    """
    Execute GSEApy annotation for a single specific cluster.

    # CORE RESPONSIBILITY
    This tool provides detailed functional enrichment analysis for ONE cluster:
    1. Extracts marker genes for the specified cluster from DGE results
    2. Queries Enrichr databases to identify cell types or pathways
    3. Returns detailed enrichment results with statistics

    # WHEN TO USE
    - After DGE analysis is complete
    - When user asks about a specific cluster (e.g., "What is cluster 0?")
    - When you need detailed annotation for one cluster
    - After batch annotation, to get more details on a specific cluster

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: DGE analysis must be completed. The `snapshot_id` should point to a DGE result.
    2. **Cluster ID**: Must be a valid cluster ID from the clustering results (e.g., "0", "1", "2").
    3. **Categories**: Choose specific databases based on analysis goal:
       - "cellmarker": For cell type identification
       - "functional": For pathway analysis (GO, KEGG)
       - "msigdb": For hallmark gene sets

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a DGE analysis snapshot.
    - **`cluster_id` is MANDATORY**: Must be a valid cluster ID (string or number).

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from DGE analysis result (e.g., 's_dge_wilcoxon_...'). REQUIRED.
        cluster_id: The specific cluster ID to annotate (e.g., "0", "1", "2"). REQUIRED.
        top_n_genes: Number of top marker genes to use for enrichment (default: 100).
        categories: List of database categories to query. Options: ["cellmarker", "functional", "msigdb"].
                   If None, uses all available databases.

    Returns:
        Detailed annotation results for the specified cluster including:
        - Predicted cell type or enriched pathways
        - Statistical significance scores
        - Top enriched terms
        - New snapshot ID
    """
    try:
        payload = {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "top_n_genes": top_n_genes,
            "cluster_id": str(cluster_id)  # Ensure it's a string
        }
        if categories:
            payload["categories"] = categories

        data = await call_backend("POST", "/annotation/gseapy/single", payload)

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        enrichment_results = data.get("enrichment_results", {})
        cluster_id_result = data.get("cluster_id", cluster_id)
        msg = data.get("msg", "")

        # Format enrichment results
        results_summary = []
        if isinstance(enrichment_results, dict):
            # Try to extract top enriched terms
            for key, value in list(enrichment_results.items())[:5]:
                if isinstance(value, dict):
                    term = value.get("Term", key)
                    pvalue = value.get("P-value", "N/A")
                    results_summary.append(f"  - {term}: p-value = {pvalue}")
                else:
                    results_summary.append(f"  - {key}: {value}")

        return (
            f"Single Cluster GSEApy Annotation Complete!\n"
            f"\n"
            f"- Cluster ID: {cluster_id_result}\n"
            f"- Top genes used: {top_n_genes}\n"
            f"- Databases: {', '.join(categories) if categories else 'All available'}\n"
            f"\n"
            f"Enrichment Results:\n"
            f"{chr(10).join(results_summary) if results_summary else '  (Processing results...)'}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Message: {msg}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review enrichment results to understand cluster identity\n"
            f"  2. Compare with other clusters if needed\n"
            f"  3. Consider running batch annotation for all clusters"
        )
    except Exception as e:
        return f"Single cluster GSEApy annotation failed: {str(e)}"


@tool(
    description=(
            "Run CellTypist automated cell type annotation for the entire dataset. "
            "Use this for machine learning-based cell type prediction using pre-trained models."
    ),
    parse_docstring=True
)
async def run_celltypist_annotation(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        model_names: List[str] = None,
        majority_voting: bool = True,
        target_cluster_col: str = "leiden"
) -> str:
    """
    Execute CellTypist automated cell type annotation using pre-trained models.

    # CORE RESPONSIBILITY
    This tool uses machine learning models to predict cell types:
    1. Loads clustering results (with cell embeddings)
    2. Applies pre-trained CellTypist models to predict cell identities
    3. Optionally performs majority voting at cluster level
    4. Assigns dominant cell type label to each cluster

    # WHEN TO USE
    - After clustering is complete
    - When user wants ML-based cell type prediction
    - As an alternative to GSEApy annotation
    - When working with immune cells (many models available)

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: Clustering must be completed. The `snapshot_id` should point to a clustering result.
    2. **Model Selection**: Choose models based on cell type:
       - "Immune_All_Low.pkl": General immune cells (default)
       - "Human_Lung_Atlas.pkl": Lung-specific
       - "Human_All.pkl": All human cell types
    3. **Majority Voting**: If True, aggregates predictions at cluster level (recommended).
    4. **Target Cluster Column**: The clustering column name (usually "leiden" or "louvain").

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a clustering snapshot (e.g., 's_cluster_leiden_...').
    - **Clustering Required**: This tool requires clustering to be completed first.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from clustering result (e.g., 's_cluster_leiden_...'). REQUIRED.
        model_names: List of CellTypist model names to apply. Default: ["Immune_All_Low.pkl"].
        majority_voting: Whether to perform cluster-level majority voting (default: True).
        target_cluster_col: Column name for clustering (e.g., "leiden", "louvain"). Default: "leiden".

    Returns:
        Summary of CellTypist annotation results including:
        - Number of cells annotated
        - Predicted cell types per cluster
        - Model confidence scores
        - New snapshot ID
    """
    try:
        if model_names is None:
            model_names = ["Immune_All_Low.pkl"]

        data = await call_backend("POST", "/annotation/celltypist", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "model_names": model_names,
            "majority_voting": majority_voting,
            "target_cluster_col": target_cluster_col
        })

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        enrichment_results = data.get("enrichment_results", {})
        msg = data.get("msg", "")

        # Format annotation summary
        summary_lines = []
        if isinstance(enrichment_results, dict):
            # Try to extract cluster-level predictions
            for key, value in list(enrichment_results.items())[:10]:
                if isinstance(value, dict):
                    cell_type = value.get("predicted_cell_type", value.get("cell_type", "Unknown"))
                    confidence = value.get("confidence", "N/A")
                    summary_lines.append(f"  Cluster {key}: {cell_type} (confidence: {confidence})")
                else:
                    summary_lines.append(f"  {key}: {value}")

        return (
            f"CellTypist Annotation Complete!\n"
            f"\n"
            f"- Models used: {', '.join(model_names)}\n"
            f"- Majority voting: {'Enabled' if majority_voting else 'Disabled'}\n"
            f"- Cluster column: {target_cluster_col}\n"
            f"\n"
            f"Cell Type Predictions:\n"
            f"{chr(10).join(summary_lines) if summary_lines else '  (Processing results...)'}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Message: {msg}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review predicted cell types for each cluster\n"
            f"  2. Compare with GSEApy annotation results if available\n"
            f"  3. Verify predictions match expected cell types"
        )
    except Exception as e:
        return f"CellTypist annotation failed: {str(e)}"
