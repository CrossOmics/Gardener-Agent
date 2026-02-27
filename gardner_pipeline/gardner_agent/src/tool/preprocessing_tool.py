from typing import Optional, Literal, Dict
from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description="Run the complete preprocessing pipeline (QC -> Filter -> HVG -> PCA -> Neighbors) in one call. "
                "Use this to transform raw genomic data into a graph-based representation for clustering.",
    parse_docstring=True
)
async def run_full_preprocessing(
        project_id: str,
        dataset_id: str,
        organism: str = "Human",
        min_genes: int = 200,
        min_cells: int = 3,
        pct_mt_max: Optional[float] = None,
        pct_hb_max: Optional[float] = None,
        n_top_genes: int = 2000,
        flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = "seurat",
        target_sum: float = 1e4,
        n_comps: int = 50,
        n_neighbors: int = 15,
        n_pcs: int = 30,
        custom_prefixes: Optional[Dict[str, str]] = None,
        skip_qc_calculation: bool = False,
        skip_qc_filter: bool = False,
        skip_hvg: bool = False,
        skip_pca: bool = False,
        skip_neighbors: bool = False
) -> str:
    """
    Execute the full single-cell preprocessing workflow from raw data.

    # CORE RESPONSIBILITY
    Use this tool when a raw dataset has just been imported and needs initial processing.
    It sequentially handles Quality Control, Normalization, Feature Selection, and Dimensionality Reduction.

    # STRATEGY & BEST PRACTICES
    1. **Starting Point**: No snapshot_id is needed; it starts from the raw dataset_id.
    2. **QC Thresholds**:
       - Human: mitochondrial genes usually start with 'MT-'.
       - Strict Filtering: set `pct_mt_max` to 5.0 or 10.0 to remove low-quality/dying cells.
    3. **HVG**: 2000 genes is standard. For complex tissues, consider increasing to 3000.
    4. **Neighbors**: n_neighbors=15 and n_pcs=30 is the robust default for subsequent Leiden clustering.

    # CONSTRAINTS
    - `project_id` and `dataset_id` are MANDATORY.
    - If user asks for "standard preprocessing", keep defaults but ensure `organism` matches.

    Args:
        project_id: Unique business ID of the project. REQUIRED.
        dataset_id: Unique ID of the raw dataset. REQUIRED.
        organism: Biological source (e.g., "Human", "Mouse").
        min_genes: Minimum genes per cell to keep (default: 200).
        min_cells: Minimum cells expressing a gene to keep (default: 3).
        pct_mt_max: Max allowed mitochondrial % per cell (e.g., 5.0).
        pct_hb_max: Max allowed hemoglobin % per cell.
        n_top_genes: Number of highly variable genes to select (default: 2000).
        flavor: Selection algorithm: 'seurat', 'cell_ranger', or 'seurat_v3'.
        target_sum: Normalization factor (default: 10000).
        n_comps: PCA components to compute (default: 50).
        n_neighbors: Neighbors for the adjacency graph (default: 15).
        n_pcs: Number of PCs used for neighbor searching (default: 30).
        custom_prefixes: Map for gene prefixes, e.g., {"mt": "MT-", "hb": "HB-"}.
        skip_qc_calculation: Bypass metric calculation.
        skip_qc_filter: Bypass filtering steps.
        skip_hvg: Bypass HVG selection.
        skip_pca: Bypass PCA.
        skip_neighbors: Bypass Graph construction.
    """
    try:
        payload = {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "organism": organism,
            "custom_prefixes": custom_prefixes,
            "min_genes": min_genes,
            "min_cells": min_cells,
            "pct_mt_max": pct_mt_max,
            "pct_hb_max": pct_hb_max,
            "n_top_genes": n_top_genes,
            "flavor": flavor,
            "target_sum": target_sum,
            "n_comps": n_comps,
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "skip_qc_calculation": skip_qc_calculation,
            "skip_qc_filter": skip_qc_filter,
            "skip_hvg": skip_hvg,
            "skip_pca": skip_pca,
            "skip_neighbors": skip_neighbors
        }

        # Remove None values to use backend defaults
        payload = {k: v for k, v in payload.items() if v is not None}

        data = await call_backend("POST", "/preprocessing/run/full", payload)

        executed = data.get("executed_steps", [])
        final_stage = data.get("final_stage", "unknown")
        # In the context of our agent system, the anndata_cache_path usually
        # aligns with the Snapshot ID for subsequent tools.
        snapshot_id = data.get("snapshot_id", "N/A")

        return (
            f"Preprocessing Pipeline Successful.\n"
            f"- Executed Steps: {', '.join(executed)}\n"
            f"- Final Stage Reached: {final_stage}\n"
            f"- Generated Snapshot ID: {snapshot_id}\n"
            f"- Data Path: {data.get('anndata_cache_path')}\n"
            f"Observation: Preprocessing is finished. You can now proceed to clustering or visualize the QC plots."
        )
    except Exception as e:
        return f"Full preprocessing pipeline failed: {str(e)}"


@tool(
    description=(
            "Calculate Quality Control (QC) metrics for raw data. "
            "Use this to assess data quality before filtering."
    ),
    parse_docstring=True
)
async def calculate_qc_metrics(
        project_id: str,
        dataset_id: str,
        organism: str = "Human"
) -> str:
    """
    Calculate Quality Control metrics for raw single-cell data.

    # CORE RESPONSIBILITY
    This tool computes QC metrics including:
    1. Total counts per cell
    2. Number of genes per cell
    3. Mitochondrial gene percentage
    4. Hemoglobin gene percentage
    5. Generates visualization plots (violin plots, scatter plots)

    # WHEN TO USE
    - Before applying QC filters
    - When user wants to assess data quality
    - When user asks: "Check data quality" or "Calculate QC metrics"
    - As the first step in preprocessing workflow

    # STRATEGY & BEST PRACTICES
    1. **Starting Point**: This tool starts from the raw dataset (no snapshot_id needed).
    2. **Organism Selection**: Choose "Human" or "Mouse" for automatic gene prefix detection.
    3. **Review Results**: Check the generated plots to determine appropriate filtering thresholds.

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`dataset_id` is MANDATORY**: Must refer to an uploaded raw dataset.
    - **No Snapshot Required**: This tool starts from raw data.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the raw dataset. REQUIRED.
        organism: Organism type for gene prefix detection (default: "Human"). Options: "Human", "Mouse".

    Returns:
        Summary of QC calculation results including:
        - Paths to QC visualization plots
        - Metrics JSON file path
    """
    try:
        data = await call_backend("POST", "/preprocessing/qc/calculate", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "organism": organism
        })

        # Extract key information
        metrics_json_path = data.get("metrics_json_path", "")
        violin_plot_path = data.get("violin_plot_path", "")
        scatter_mt_path = data.get("scatter_mt_path", "")
        scatter_genes_path = data.get("scatter_genes_path", "")

        return (
            f"QC Metrics Calculation Complete!\n"
            f"\n"
            f"Results:\n"
            f"  - Metrics JSON: {metrics_json_path}\n"
            f"  - Violin plot: {violin_plot_path}\n"
            f"  - Scatter plot (MT): {scatter_mt_path}\n"
            f"  - Scatter plot (Genes): {scatter_genes_path}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review QC plots to determine filtering thresholds\n"
            f"  2. Apply QC filtering with appropriate parameters\n"
            f"  3. Consider filtering based on min_genes and pct_mt_max"
        )
    except Exception as e:
        return f"QC metrics calculation failed: {str(e)}"


@tool(
    description=(
            "Apply Quality Control filtering to remove low-quality cells and genes. "
            "Use this after QC calculation to filter the dataset based on quality thresholds."
    ),
    parse_docstring=True
)
async def apply_qc_filter(
        project_id: str,
        dataset_id: str,
        min_genes: int = 200,
        min_cells: int = 3,
        pct_mt_max: Optional[float] = None,
        max_counts: Optional[int] = None
) -> str:
    """
    Apply Quality Control filtering to remove low-quality cells and genes.

    # CORE RESPONSIBILITY
    This tool filters the dataset by:
    1. Removing cells with too few genes (below min_genes threshold)
    2. Removing genes expressed in too few cells (below min_cells threshold)
    3. Optionally filtering by mitochondrial percentage (pct_mt_max)
    4. Optionally filtering by total counts (max_counts)
    5. Creates a new snapshot with filtered data

    # WHEN TO USE
    - After QC calculation is complete
    - When user wants to remove low-quality cells/genes
    - When user asks: "Filter the data" or "Apply QC filters"
    - Before HVG selection (filtering is typically done first)

    # STRATEGY & BEST PRACTICES
    1. **Starting Point**: This tool starts from the raw dataset (no snapshot_id needed).
    2. **Threshold Selection**: Based on QC plots:
       - `min_genes`: Typical range 100-500 (default 200)
       - `pct_mt_max`: Typical range 5-20% (lower = stricter)
    3. **Conservative Filtering**: Start with lenient thresholds, then tighten if needed.

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`dataset_id` is MANDATORY**: Must refer to an uploaded raw dataset.
    - **No Snapshot Required**: This tool starts from raw data.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the raw dataset. REQUIRED.
        min_genes: Minimum number of genes per cell (default: 200).
        min_cells: Minimum number of cells expressing a gene (default: 3).
        pct_mt_max: Maximum mitochondrial percentage (0-100). If None, no filtering by MT%.
        max_counts: Maximum total counts per cell. If None, no filtering by counts.

    Returns:
        Summary of filtering results including:
        - Number of cells remaining
        - Number of genes remaining
        - New snapshot ID
    """
    try:
        payload = {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "min_genes": min_genes,
            "min_cells": min_cells
        }
        if pct_mt_max is not None:
            payload["pct_mt_max"] = pct_mt_max
        if max_counts is not None:
            payload["max_counts"] = max_counts

        data = await call_backend("POST", "/preprocessing/qc/filter", payload)

        # Extract key information
        snapshot_id = data.get("snapshot_id")
        n_obs_remaining = data.get("n_obs_remaining", 0)
        n_vars_remaining = data.get("n_vars_remaining", 0)

        return (
            f"QC Filtering Complete!\n"
            f"\n"
            f"- Cells remaining: {n_obs_remaining}\n"
            f"- Genes remaining: {n_vars_remaining}\n"
            f"- Min genes per cell: {min_genes}\n"
            f"- Min cells per gene: {min_cells}\n"
            f"{f'- Max MT%: {pct_mt_max}' if pct_mt_max else ''}\n"
            f"{f'- Max counts: {max_counts}' if max_counts else ''}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Proceed to HVG selection\n"
            f"  2. Review filtering results and adjust thresholds if needed"
        )
    except Exception as e:
        return f"QC filtering failed: {str(e)}"


@tool(
    description=(
            "Select Highly Variable Genes (HVG) and perform normalization. "
            "Use this after QC filtering to identify genes with high variability for downstream analysis."
    ),
    parse_docstring=True
)
async def apply_hvg_selection(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        n_top_genes: int = 2000,
        flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = "seurat",
        target_sum: float = 1e4
) -> str:
    """
    Select Highly Variable Genes (HVG) and perform normalization.

    # CORE RESPONSIBILITY
    This tool performs:
    1. Normalization (total counts normalization)
    2. Log transformation (log1p)
    3. Highly Variable Gene identification
    4. Scaling (standardization)
    5. Creates a new snapshot with processed data

    # WHEN TO USE
    - After QC filtering is complete
    - When user wants to select variable genes for analysis
    - When user asks: "Select HVG" or "Run feature selection"
    - Before PCA (HVG selection is a prerequisite)

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: QC filtering should be completed. The `snapshot_id` should point to a filtered dataset snapshot.
    2. **Number of HVGs**: Default 2000 is usually sufficient. Adjust based on:
       - Dataset size: Larger datasets may need more HVGs
       - Analysis goal: More HVGs = more genes for analysis
    3. **Flavor Selection**: 
       - "seurat": Default, most commonly used
       - "cell_ranger": For 10x Genomics data
       - "seurat_v3": Alternative method

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a filtered dataset snapshot (e.g., 's_filter_...').
    - **QC Filtering Required**: This tool requires QC filtering to be completed first.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from QC filtering result (e.g., 's_filter_...'). REQUIRED.
        n_top_genes: Number of highly variable genes to select (default: 2000, range: 500-10000).
        flavor: HVG selection method. Options: "seurat" (default), "cell_ranger", "seurat_v3".
        target_sum: Target total counts per cell for normalization (default: 10000).

    Returns:
        Summary of HVG selection results including:
        - Number of HVGs found
        - HVG dispersion plot path
        - New snapshot ID
    """
    try:
        data = await call_backend("POST", "/preprocessing/hvg", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "n_top_genes": n_top_genes,
            "flavor": flavor,
            "target_sum": target_sum
        })

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        n_genes_found = data.get("n_genes_found", 0)
        hvg_plot_path = data.get("hvg_plot_path", "")
        msg = data.get("msg", "")

        return (
            f"HVG Selection Complete!\n"
            f"\n"
            f"- HVGs selected: {n_genes_found}\n"
            f"- Method: {flavor}\n"
            f"- Target sum: {target_sum}\n"
            f"\n"
            f"Results:\n"
            f"  - HVG dispersion plot: {hvg_plot_path}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Message: {msg}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Proceed to PCA analysis\n"
            f"  2. Review HVG dispersion plot if needed"
        )
    except Exception as e:
        return f"HVG selection failed: {str(e)}"
