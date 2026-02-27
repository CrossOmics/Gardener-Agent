from pydantic import BaseModel, Field
from typing import Optional, List, Dict

from dto.request.base_request import BaseAnalysisRequest


class RunFullAnnotationRequest(BaseAnalysisRequest):
    # GSEApy Parameters
    # If categories is None or empty, GSEApy step is skipped.
    top_n_genes: int = 100
    categories: Optional[List[str]] = Field(
        default=[],
        description="List of Enrichr libraries (e.g. ['CellMarker_2024', '']). Pass empty list or null to skip."
    )

    # CellTypist Parameters
    # If model_names is None or empty, CellTypist step is skipped.
    model_names: Optional[List[str]] = Field(
        default=[],
        description="List of CellTypist models. Pass empty list or null to skip."
    )
    majority_voting: bool = True
    target_cluster_col: Optional[str] = Field(default="leiden")


class UpdateAnnotationLabelRequest(BaseModel):
    """
    Request payload for updating specific cluster labels in an annotation snapshot.
    """
    snapshot_id: str = Field(..., description="The unique business ID of the analysis snapshot.")

    file_name: str = Field(
        ...,
        description="The key name in thumbnail_json (usually the file name without extension). e.g., 'gseapy_CellMarker_2024'"
    )

    updated_annotation: Dict[str, str] = Field(
        ...,
        description="Map of Cluster ID to New Label Name. e.g., {'0': 'New T-Cell Label', '1': 'B-Cell'}"
    )


class RunAnnotationRequest(BaseAnalysisRequest):
    top_n_genes: int = 100
    cluster_id: Optional[str] = "ALL"
    categories: Optional[List[str]] = None  # e.g. ["cellmarker", "functional"]


class RunCellTypistRequest(BaseAnalysisRequest):
    """
    Request DTO for CellTypist automated cell type annotation.
    """
    # CellTypist Parameters
    # accepts a list of model names
    model_names: List[str] = Field(
        default=["Immune_All_Low.pkl"],
        description="List of CellTypist models to apply (e.g. ['Immune_All_Low.pkl', 'Human_Lung_Atlas.pkl'])"
    )

    # Whether to refine predictions using cluster-based majority voting
    majority_voting: bool = Field(default=True)

    # Column name in adata.obs for majority voting (e.g., leiden)
    target_cluster_col: Optional[str] = Field(default="leiden")
