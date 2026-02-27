from typing import Literal

from langchain_core.tools import tool
from base.base_call import call_backend


@tool(
    description=(
        "Run Principal Component Analysis (PCA) for dimensionality reduction. "
        "Use this after HVG selection to reduce data dimensions and prepare for clustering."
    ),
    parse_docstring=True
)
async def run_pca_analysis(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        n_comps: int = 50,
        svd_solver: Literal['arpack', 'randomized', 'auto'] = "arpack"
) -> str:
    """
    Execute Principal Component Analysis (PCA) for dimensionality reduction.

    # CORE RESPONSIBILITY
    This tool performs PCA on highly variable genes to:
    1. Reduce data dimensionality (from thousands of genes to ~50 principal components)
    2. Extract variance information for Elbow Plot visualization
    3. Prepare data for neighborhood graph construction and clustering

    # WHEN TO USE
    - After HVG (Highly Variable Genes) selection is complete
    - When user wants to reduce dimensionality before clustering
    - When user asks: "Run PCA" or "Perform dimensionality reduction"
    - Before building neighborhood graph (PCA is a prerequisite)

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: HVG selection must be completed. The `snapshot_id` should point to an HVG result snapshot.
    2. **Number of Components**: Default 50 is usually sufficient. Adjust based on:
       - Dataset size: Larger datasets may need more components
       - Elbow Plot: Check variance explained to determine optimal number
    3. **SVD Solver**: 
       - "arpack": Default, good for most cases
       - "randomized": Faster for large datasets
       - "auto": Automatically selects based on data size

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be an HVG snapshot (e.g., 's_hvg_...').
    - **HVG Required**: This tool requires HVG selection to be completed first.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from HVG result (e.g., 's_hvg_...'). REQUIRED.
        n_comps: Number of principal components to compute (default: 50). Typical range: 30-100.
        svd_solver: SVD solver algorithm. Options: "arpack" (default), "randomized" (faster for large data), "auto".

    Returns:
        Summary of PCA analysis results including:
        - Variance explained by each component (for Elbow Plot)
        - PCA scatter plot path
        - New snapshot ID
    """
    try:
        data = await call_backend("POST", "/dim-reduction/pca", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "n_comps": n_comps,
            "svd_solver": svd_solver
        })

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        variance_plot_data = data.get("variance_plot_data", [])
        pca_scatter_path = data.get("pca_scatter_path", "")

        # Calculate variance explained by top components
        top_10_variance = variance_plot_data[:10] if len(variance_plot_data) >= 10 else variance_plot_data
        cumulative_variance = sum(variance_plot_data[:n_comps]) if len(variance_plot_data) >= n_comps else sum(variance_plot_data)

        return (
            f"PCA Analysis Complete!\n"
            f"\n"
            f"- Number of components: {n_comps}\n"
            f"- SVD solver: {svd_solver}\n"
            f"- Total variance explained: {cumulative_variance:.2%}\n"
            f"- Top 10 component variances: {[f'{v:.4f}' for v in top_10_variance]}\n"
            f"\n"
            f"Results:\n"
            f"  - PCA scatter plot: {pca_scatter_path}\n"
            f"  - Variance plot data: Available for Elbow Plot visualization\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Review Elbow Plot to determine optimal number of PCs\n"
            f"  2. Build neighborhood graph using selected number of PCs\n"
            f"  3. Proceed to clustering analysis"
        )
    except Exception as e:
        return f"PCA analysis failed: {str(e)}"


@tool(
    description=(
        "Build neighborhood graph for cell-cell relationships. "
        "Use this after PCA to construct the k-nearest neighbor graph required for UMAP and clustering."
    ),
    parse_docstring=True
)
async def build_neighborhood_graph(
        project_id: str,
        dataset_id: str,
        snapshot_id: str,
        n_neighbors: int = 15,
        n_pcs: int = 30
) -> str:
    """
    Construct k-nearest neighbor graph for downstream analysis.

    # CORE RESPONSIBILITY
    This tool builds the cell-cell neighbor graph by:
    1. Computing k-nearest neighbors based on PCA space
    2. Creating connectivity matrix for UMAP embedding
    3. Preparing graph structure for clustering algorithms (Leiden/Louvain)

    # WHEN TO USE
    - After PCA analysis is complete
    - When user wants to prepare data for clustering or UMAP
    - When user asks: "Build neighbors graph" or "Prepare for clustering"
    - Before running clustering (neighborhood graph is required)

    # STRATEGY & BEST PRACTICES
    1. **Prerequisites**: PCA must be completed. The `snapshot_id` should point to a PCA result snapshot.
    2. **Number of Neighbors**: Default 15 is usually good. Adjust based on:
       - Higher (20-30): Smoother clusters, fewer small groups
       - Lower (10-15): More fine-grained clusters, may fragment
    3. **Number of PCs**: Based on Elbow Plot results:
       - Use PCs that explain most variance (typically 20-50)
       - Too few PCs: May lose important structure
       - Too many PCs: May include noise

    # CONSTRAINTS
    - **`project_id` is MANDATORY**: Use the exact ID from conversation context. Never guess.
    - **`snapshot_id` is MANDATORY**: Must be a PCA snapshot (e.g., 's_pca_...').
    - **PCA Required**: This tool requires PCA analysis to be completed first.

    Args:
        project_id: Unique business ID of the project (e.g., 'p_2023...'). REQUIRED.
        dataset_id: Unique ID of the dataset.
        snapshot_id: Snapshot ID from PCA result (e.g., 's_pca_...'). REQUIRED.
        n_neighbors: Number of nearest neighbors to connect (default: 15). Typical range: 10-30.
        n_pcs: Number of principal components to use for neighbor search (default: 30). Based on Elbow Plot.

    Returns:
        Summary of neighborhood graph construction including:
        - Number of connections created
        - New snapshot ID
    """
    try:
        data = await call_backend("POST", "/dim-reduction/neighbors", {
            "project_id": project_id,
            "dataset_id": dataset_id,
            "snapshot_id": snapshot_id,
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs
        })

        # Extract key information
        snapshot_id_new = data.get("snapshot_id")
        n_connectivities = data.get("n_connectivities", 0)

        return (
            f"Neighborhood Graph Construction Complete!\n"
            f"\n"
            f"- Number of neighbors: {n_neighbors}\n"
            f"- Principal components used: {n_pcs}\n"
            f"- Total connections: {n_connectivities}\n"
            f"\n"
            f"New Snapshot ID: {snapshot_id_new}\n"
            f"\n"
            f"Next Steps:\n"
            f"  1. Proceed to clustering analysis (Leiden/Louvain)\n"
            f"  2. The graph is ready for UMAP visualization\n"
            f"  3. Adjust n_neighbors if clustering results are unsatisfactory"
        )
    except Exception as e:
        return f"Neighborhood graph construction failed: {str(e)}"