from langchain_core.tools import tool
from typing import Optional, List, Dict, Any
from base.base_call import call_backend


@tool(
    description=(
            "Execute the complete end-to-end single-cell RNA-seq analysis pipeline. "
            "Runs: Preprocessing -> Clustering -> DEG -> Annotation in one go."
    ),
    parse_docstring=True
)
async def run_whole_pipeline(
        project_id: str,
        dataset_id: str,
        # Clustering
        clustering_method: str,
        # Preprocessing
        organism: str = "Human",
        min_genes: int = 200,
        max_genes: Optional[int] = None,
        min_cells: int = 3,
        pct_mt_max: Optional[float] = 5.0,
        pct_hb_max: Optional[float] = None,
        n_top_genes: int = 2000,
        flavor: str = "seurat",
        target_sum: float = 1e4,
        n_comps: int = 50,
        n_neighbors: int = 15,
        n_pcs: int = 30,

        clustering_resolution: float = 0.5,
        # DEG
        deg_method: str = "wilcoxon",
        deg_n_top_genes: int = 25,
        # Annotation
        annotation_categories: Optional[List[str]] = None,
        celltypist_model: Optional[str] = "Immune_All_Low.pkl",
        annotation_top_n_genes: int = 100,
        majority_voting: bool = True
) -> str:
    """
    Execute the complete single-cell analysis pipeline from raw data to annotated clusters.

    # CORE RESPONSIBILITY
    Use this tool when the user wants to "run everything" or "analyze this dataset" without
    specifying individual steps. It automates:
    1. Preprocessing (QC, Filtering, Normalization, PCA, Neighbors)
    2. Clustering (leiden/louvain/cplearn)
    3. Differential Expression (Marker Genes)
    4. Cell Type Annotation (GSEApy & CellTypist)

    # STRATEGY
    - Starts from the RAW dataset (no snapshot_id needed).
    - Returns a final snapshot ID that is fully annotated and ready for visualization.

    Args:
        project_id: Unique project identifier. REQUIRED.
        dataset_id: Dataset to analyze. REQUIRED.
        organism: "Human" or "Mouse" (used for QC and Annotation).
        min_genes: QC filter: min genes per cell (default: 200).
        max_genes: QC filter: max genes per cell (default: None).
        min_cells: QC filter: min cells per gene (default: 3).
        pct_mt_max: QC filter: max mitochondrial % (default: 5.0).
        pct_hb_max: QC filter: max hemoglobin % (default: None).
        n_top_genes: HVG selection count (default: 2000).
        flavor: HVG selection method ('seurat', 'cell_ranger', 'seurat_v3').
        target_sum: Target sum for normalization (default: 1e4).
        n_comps: Number of PCA components (default: 50).
        n_neighbors: Number of neighbors for graph (default: 15).
        n_pcs: Number of PCs to use for neighbors (default: 30).
        clustering_method: Clustering algorithm ('leiden', 'louvain', 'cplearn').
        clustering_resolution: Resolution for clustering (default: 0.5).
        deg_method: Test method for marker genes ('wilcoxon', 't-test', 'logreg').
        deg_n_top_genes: Number of top marker genes to identify per cluster (default: 25).
        annotation_categories: List of Enrichr libraries (default: ["CellMarker_2024"]).
        celltypist_model: Model for automated annotation (default: "Immune_All_Low.pkl"). Pass None to skip.
        annotation_top_n_genes: Number of genes to use for GSEApy annotation (default: 100).
        majority_voting: Whether to refine CellTypist predictions with majority voting (default: True).
    """
    try:
        # 1. Construct the nested payload matching RunWholePipelineRequest
        payload = {
            "project_id": project_id,
            "dataset_id": dataset_id,

            # Nested Preprocessing DTO
            "preprocessing_params": {
                "organism": organism,
                "min_genes": min_genes,
                "max_genes": max_genes,
                "min_cells": min_cells,
                "pct_mt_max": pct_mt_max,
                "pct_hb_max": pct_hb_max,
                "n_top_genes": n_top_genes,
                "flavor": flavor,
                "target_sum": target_sum,
                "n_comps": n_comps,
                "n_neighbors": n_neighbors,
                "n_pcs": n_pcs,
                # Defaulting others
                "svd_solver": "arpack"
            },

            # Nested Clustering DTO
            "clustering_params": {
                "method": clustering_method,
                "resolution": clustering_resolution,
                "run_hierarchical": True
            },

            # Nested DEG DTO
            "deg_params": {
                "method": deg_method,
                "groupby": clustering_method,
                "n_top_genes": deg_n_top_genes,
                "use_raw": True
            },

            # Nested Annotation DTO
            "annotation_params": {
                "categories": annotation_categories or ["CellMarker_2024"],
                "model_names": [celltypist_model] if celltypist_model else [],
                "top_n_genes": annotation_top_n_genes,
                "majority_voting": majority_voting,
                "target_cluster_col": clustering_method
            }
        }

        # 2. Call Backend
        data = await call_backend("POST", "/pipeline/run", payload)

        # 3. Parse Response (RunWholePipelineResponse)
        pipeline_id = data.get("pipeline_run_id")
        final_snapshot = data.get("final_snapshot_id")
        status_msg = data.get("msg", "Completed")
        snapshots_map = data.get("snapshots", {})

        # Format a readable summary of intermediate snapshots
        steps_summary = "\n".join([f"- {k}: {v}" for k, v in snapshots_map.items()])

        return (
            f"Whole Pipeline Execution Successful!\n"
            f"Run ID: {pipeline_id}\n"
            f"Status: {status_msg}\n"
            f"Final Snapshot ID: {final_snapshot}\n"
            f"\n"
            f"Generated Snapshots:\n"
            f"{steps_summary}\n"
            f"\n"
            f"Observation: The dataset has been preprocessed, clustered, and annotated. "
            f"You can now visualize the results using the Final Snapshot ID."
        )

    except Exception as e:
        return f"Pipeline execution failed: {str(e)}"
