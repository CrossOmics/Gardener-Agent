from pydantic import Field
from typing import Optional, Literal, Dict

from dto.request.base_request import BaseAnalysisRequest


class FullPreprocessingRequest(BaseAnalysisRequest):
    """
    Request DTO for running the full preprocessing pipeline.

    Pipeline:
    QC Calculation -> QC Filter -> HVG -> PCA -> Neighbors

    The entire pipeline operates on a single AnnData object.
    Intermediate steps do NOT persist snapshots.
    """

    # QC Calculation Parameters
    organism: Optional[str] = Field("Human", description="Organism type (Human/Mouse) for gene prefix detection")
    custom_prefixes: Optional[Dict[str, str]] = Field(None, description="Custom gene prefixes (e.g., {'mt': 'MT-'})")

    # Skip switches (all optional, default False)
    skip_qc_calculation: bool = Field(False, description="Skip QC metrics calculation")
    skip_qc_filter: bool = Field(False, description="Skip QC filtering")
    skip_hvg: bool = Field(False, description="Skip HVG selection")
    skip_pca: bool = Field(False, description="Skip PCA computation")
    skip_neighbors: bool = Field(False, description="Skip neighbor graph construction")

    # QC Filter — Cell-level (filter_cells)
    min_genes: int = Field(200, description="Minimum number of genes per cell")
    cell_min_counts: Optional[int] = Field(None, description="Minimum total counts per cell")
    max_genes: Optional[int] = Field(None, description="Maximum number of genes per cell")
    cell_max_counts: Optional[int] = Field(None, description="Maximum total counts per cell")

    # QC Filter — Gene-level (filter_genes)
    min_cells: int = Field(3, description="Minimum number of cells expressing a gene")
    gene_min_counts: Optional[int] = Field(None, description="Minimum total counts per gene")
    max_cells: Optional[int] = Field(None, description="Maximum number of cells expressing a gene")
    gene_max_counts: Optional[int] = Field(None, description="Maximum total counts per gene")

    # QC Metric-based Filters (require QC Calculation)
    pct_mt_max: Optional[float] = Field(None, description="Maximum mitochondrial percentage (0–100)")
    pct_hb_max: Optional[float] = Field(None, description="Maximum hemoglobin percentage (0–100)")

    # HVG Parameters
    n_top_genes: int = Field(2000, description="Number of HVGs to select")
    flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = Field('seurat', description="HVG selection method")
    target_sum: float = Field(1e4, description="Target total counts per cell for normalization")

    # PCA Parameters
    n_comps: int = Field(50, description="Number of PCA components")
    svd_solver: str = Field('arpack', description="SVD solver for PCA")

    # Neighbors Parameters
    n_neighbors: int = Field(15, description="Number of neighbors")
    n_pcs: int = Field(30, description="Number of PCs used to build neighbor graph")
