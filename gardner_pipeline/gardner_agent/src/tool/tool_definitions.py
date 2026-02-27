from tool.clustering_tool import run_clustering, merge_clusters, run_sub_clustering
from tool.whole_pipeline_tool import run_whole_pipeline
from tool.dge_tool import run_dge_analysis
from tool.annotation_tool import run_full_annotation_pipeline, update_annotation_labels
from tool.dim_reduction_tool import run_pca_analysis, build_neighborhood_graph
from tool.preprocessing_tool import run_full_preprocessing, calculate_qc_metrics, apply_qc_filter, apply_hvg_selection
from tool.snapshot_tool import get_snapshot_details, delete_snapshot, query_latest_stage_snapshot
from tool.project_management_tool import *
from tool.user_preference_tool import *
from tool.search_annotation_method_tool import *

# Mapping Registry: Tool Name -> Function Object
TOOL_REGISTRY = {
    "run_clustering": run_clustering,
    "run_whole_pipeline": run_whole_pipeline,
    "merge_clusters": merge_clusters,
    "run_sub_clustering": run_sub_clustering,
    "run_dge_analysis": run_dge_analysis,
    "run_pca_analysis": run_pca_analysis,
    "build_neighborhood_graph": build_neighborhood_graph,
    "run_full_preprocessing": run_full_preprocessing,
    "calculate_qc_metrics": calculate_qc_metrics,
    "apply_qc_filter": apply_qc_filter,
    "apply_hvg_selection": apply_hvg_selection,
    "delete_snapshot": delete_snapshot,
    "query_latest_stage_snapshot": query_latest_stage_snapshot,
    "get_snapshot_details": get_snapshot_details,
    "create_new_project": create_new_project,
    "get_workspace_root": get_workspace_root,
    "delete_dataset": delete_dataset,
    "run_full_annotation_pipeline": run_full_annotation_pipeline,
    "update_annotation_labels": update_annotation_labels,
    "create_user_preference": create_user_preference,
    "get_user_preference": get_user_preference,
    "get_latest_preference": get_latest_preference,
    "update_user_preference": update_user_preference,
    "search_annotation_methods": search_annotation_methods,
}

"""
Tool Metadata Registry
This dictionary serves as the 'Menu' (Level 1) for the Planner/Router LLM.
It maps the function name to a high-level description of when to use it.
"""
"""
Tool Metadata Registry
This dictionary serves as the 'Menu' (Level 1) for the Planner/Router LLM.
It maps the function name to a highly detailed specification of its purpose, usage, and outputs.
"""

TOOLS_METADATA = {
    # ==========================================
    # 1. End-to-End & Macro Pipelines
    # ==========================================
    "run_whole_pipeline": (
        "[Type]: Execution\n"
        "[Description]: Executes the complete end-to-end single-cell RNA-seq analysis pipeline from raw data to final cell type annotation.\n"
        "[When to use]: When the user says 'run the standard pipeline', 'process this dataset from scratch', or wants to go from raw counts to annotated clusters in one go.\n"
        "[How to use (Hint)]: Pass default or user-specified thresholds (e.g., min_genes, target_sum, resolution, annotation models) in the suggestion.\n"
        "[Expected Output]: A final Annotation Snapshot ID, a comprehensive JSON summary report, and multiple QC/UMAP/Annotation plots."
    ),
    "run_full_preprocessing": (
        "[Type]: Execution\n"
        "[Description]: Runs the complete preprocessing workflow (QC Calculation -> QC Filter -> HVG -> PCA -> Neighborhood Graph) in a single sequence.\n"
        "[When to use]: When the user wants to prepare raw data for clustering but hasn't specified step-by-step instructions. (e.g., 'Preprocess the data', 'Run QC and PCA').\n"
        "[How to use (Hint)]: Extract any mentioned QC thresholds (pct_mt, min_cells) and PCA parameters from user input for the suggestion.\n"
        "[Expected Output]: A Preprocessing Snapshot ID, ready to be fed into clustering tools, along with variance/scatter plots."
    ),

    # ==========================================
    # 2. Step-by-Step Preprocessing
    # ==========================================
    "calculate_qc_metrics": (
        "[Type]: Execution (Non-destructive)\n"
        "[Description]: Calculates Quality Control (QC) metrics (e.g., mitochondrial percentage, gene counts) without removing any cells.\n"
        "[When to use]: When the user wants to assess data health BEFORE deciding on filtering thresholds. (e.g., 'Check the QC metrics', 'Show me the violin plots for raw data').\n"
        "[How to use (Hint)]: Needs Dataset ID and current Snapshot ID.\n"
        "[Expected Output]: A new Snapshot ID containing calculated QC columns in .obs, and QC violin plots."
    ),
    "apply_qc_filter": (
        "[Type]: Execution\n"
        "[Description]: Hard-filters the dataset by removing low-quality cells and uninformative genes based on specific thresholds.\n"
        "[When to use]: When the user explicitly sets limits (e.g., 'Filter out cells with > 5% MT', 'Remove cells with less than 200 genes').\n"
        "[How to use (Hint)]: Specify exact thresholds in the suggestion (e.g., `pct_mt_max`: 5.0, `min_genes`: 200).\n"
        "[Expected Output]: A filtered Snapshot ID and a summary report of how many cells/genes were retained vs. removed."
    ),
    "apply_hvg_selection": (
        "[Type]: Execution\n"
        "[Description]: Normalizes the data and selects Highly Variable Genes (HVG).\n"
        "[When to use]: When the user asks to 'normalize data', 'find variable genes', or 'prepare for PCA'.\n"
        "[How to use (Hint)]: Suggest target_sum (default 10000) and n_top_genes (default 2000).\n"
        "[Expected Output]: A Snapshot ID containing normalized data and HVG dispersion plots."
    ),
    "run_pca_analysis": (
        "[Type]: Execution\n"
        "[Description]: Performs Principal Component Analysis (PCA) for dimensionality reduction.\n"
        "[When to use]: When the user says 'Run PCA' or 'Calculate top principal components'.\n"
        "[How to use (Hint)]: Extract `n_pcs` (default 30-50) if mentioned by the user.\n"
        "[Expected Output]: A Snapshot ID with PCA embeddings and an elbow/variance plot."
    ),
    "build_neighborhood_graph": (
        "[Type]: Execution\n"
        "[Description]: Computes the k-nearest neighbor (KNN) graph and generates the UMAP embedding based on PCA results.\n"
        "[When to use]: When the user asks to 'build the graph', 'generate UMAP', or 'prepare for leiden(or louvain / cplearn) clustering'.\n"
        "[How to use (Hint)]: Needs `n_neighbors` and `n_pcs` to use.\n"
        "[Expected Output]: A Snapshot ID containing the connectivities graph and UMAP coordinates."
    ),

    # ==========================================
    # 3. Clustering & Modification
    # ==========================================
    "run_clustering": (
        "[Type]: Execution\n"
        "[Description]: Executes community detection clustering (leiden(or louvain / cplearn)) on the neighborhood graph.\n"
        "[When to use]: When the user asks to 'cluster the data', 'run leiden(or louvain / cplearn)', or wants to change the global resolution (e.g., 'increase resolution to 1.2 to split clusters').\n"
        "[How to use (Hint)]: Explicitly define the `resolution` parameter in the suggestion based on user intent (higher = more clusters, lower = fewer clusters).\n"
        "[Expected Output]: A Clustering Snapshot ID with a new UMAP plot colored by clusters."
    ),
    "run_sub_clustering": (
        "[Type]: Execution\n"
        "[Description]: Isolates specific clusters and re-clusters ONLY those cells at a higher resolution to reveal finer substructures.\n"
        "[When to use]: When a user identifies a cluster is too broad or mixed (e.g., 'Sub-cluster cluster 0', 'Break down the T cell cluster').\n"
        "[How to use (Hint)]: Must provide the target cluster ID(s) and a local resolution in the suggestion.\n"
        "[Expected Output]: A new Snapshot ID where ONLY the target cluster is split into sub-clusters (e.g., 0_0, 0_1), leaving others intact."
    ),
    "merge_clusters": (
        "[Type]: Execution\n"
        "[Description]: Manually forces two or more existing clusters to merge into a single group.\n"
        "[When to use]: When the user spots over-clustering and wants precise, localized corrections without changing the global resolution (e.g., 'Merge cluster 2 and 5').\n"
        "[How to use (Hint)]: Extract the exact integer IDs of the clusters to merge and provide them as a list in the suggestion.\n"
        "[Expected Output]: A new Snapshot ID, updated cluster metadata, and a refreshed UMAP plot reflecting the merged group."
    ),

    # ==========================================
    # 4. Differential Expression & Annotation
    # ==========================================
    "run_dge_analysis": (
        "[Type]: Execution\n"
        "[Description]: Calculates Differential Gene Expression (DGE) using statistical tests (e.g., Wilcoxon) to find marker genes defining each cluster.\n"
        "[When to use]: When the user asks 'What are the marker genes for cluster 1?', 'Find DEGs', or 'What defines these clusters?'.\n"
        "[How to use (Hint)]: Suggest the statistical method (`wilcoxon` or `t-test`) and how many top genes to return.\n"
        "[Expected Output]: A Snapshot ID containing rank_genes_groups data and marker gene dot plots/heatmaps."
    ),
    "run_full_annotation_pipeline": (
        "[Type]: Execution\n"
        "[Description]: Automatically assigns biological cell type identities to clusters using machine learning (CellTypist) and enrichment (GSEApy).\n"
        "[When to use]: When the user asks 'Annotate the clusters', 'What cell types are present?', or 'Run CellTypist'.\n"
        "[How to use (Hint)]: If the user mentions a specific tissue or model (e.g., 'Immune system'), suggest the corresponding `model_names`.\n"
        "[Expected Output]: An Annotation Snapshot ID, comprehensive cell type prediction plots, and a JSON confidence report."
    ),
    "update_annotation_labels": (
        "[Type]: Execution\n"
        "[Description]: Manually overrides or refines the biological label of specific clusters.\n"
        "[When to use]: When the user wants to correct the AI annotation (e.g., 'Rename cluster 3 to B-Cells', 'Cluster 0 should be Macrophages').\n"
        "[How to use (Hint)]: Provide a mapping dictionary in the suggestion (e.g., {'3': 'B-Cells'}).\n"
        "[Expected Output]: A Snapshot ID with updated `.obs` labels and a refreshed UMAP."
    ),

    # ==========================================
    # 5. Context Retrieval & Project Management
    # ==========================================
    "query_latest_stage_snapshot": (
        "[Type]: Context Retrieval\n"
        "[Description]: Looks up the most recent snapshot ID for a given analysis stage (e.g., 'Preprocessing', 'Clustering').\n"
        "[When to use]: IMPORTANT: Use this ONLY if the user asks to perform an action but the system context lacks a valid `snapshot_id`. Helps locate the data to resume work.\n"
        "[How to use (Hint)]: Pass the target stage name.\n"
        "[Expected Output]: Returns a Snapshot ID string directly to the Agent memory (No new snapshot created)."
    ),
    "get_snapshot_details": (
        "[Type]: Context Retrieval\n"
        "[Description]: Retrieves the hyperparameters, metadata, and cluster stats of a specific snapshot.\n"
        "[When to use]: When the user asks 'What parameters did we use for this result?', 'How many cells are in cluster 1?', or to establish the baseline before re-running an analysis.\n"
        "[Expected Output]: A JSON dictionary of metadata (No new snapshot created)."
    ),
    "search_annotation_methods": (
        "[Type]: Context Retrieval\n"
        "[Description]: Searches the system for available GSEApy databases or CellTypist pre-trained models.\n"
        "[When to use]: When the user asks 'What annotation models are available?', 'Can you annotate brain tissue?'.\n"
        "[Expected Output]: A text list of valid model names to be used in the annotation pipeline."
    ),
    "create_new_project": (
        "[Type]: Context/Execution\n"
        "[Description]: Initializes a fresh project workspace.\n"
        "[When to use]: 'Create a new project named X'.\n"
        "[Expected Output]: A new Project ID."
    ),
    "delete_dataset": (
        "[Type]: Execution (Destructive)\n"
        "[Description]: Deletes a dataset and associated files from the disk/database.\n"
        "[When to use]: 'Delete this dataset'. Can accept a `keep_final=True` parameter to prune intermediate files and keep only the final result.\n"
        "[Expected Output]: Deletion confirmation message."
    ),
    "delete_snapshot": (
        "[Type]: Execution (Destructive)\n"
        "[Description]: Permanently removes a specific snapshot and its physical `.h5ad` file.\n"
        "[When to use]: 'Delete the failed clustering attempt', 'Remove this snapshot'.\n"
        "[Expected Output]: Deletion confirmation message."
    ),

    # ==========================================
    # 6. User Preferences & Settings
    # ==========================================
    "create_user_preference": "[Type]: Settings\n[Description]: Save current pipeline parameter thresholds as a named profile.",
    "get_user_preference": "[Type]: Settings\n[Description]: Load a previously saved configuration profile.",
    "get_latest_preference": "[Type]: Settings\n[Description]: Retrieve the last used configuration parameters.",
    "update_user_preference": "[Type]: Settings\n[Description]: Modify an existing parameter profile.",
    "get_workspace_root": "[Type]: Context\n[Description]: Retrieve the absolute local path for the data storage root."
}


def get_tool_definitions_str() -> str:
    """
    Helper to format the dictionary into a string for the LLM System Prompt.
    Format:
    - tool_name: description
    """
    lines = []
    for name, desc in TOOLS_METADATA.items():
        lines.append(f"- {name}: {desc}")
    return "\n".join(lines)
