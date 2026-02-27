from typing import List
from pydantic import Field

from dto.request.base_request import BaseAnalysisRequest


class SubClusteringRequest(BaseAnalysisRequest):
    """
    Request DTO for running sub-clustering on specific clusters.
    """
    target_clusters: List[str] = Field(..., description="List of cluster IDs to re-cluster (e.g., ['1', '3'])")

    method: str = Field("leiden", description="Clustering algorithm (leiden, louvain)")
    resolution: float = Field(0.5, description="Resolution for the sub-clustering")

    overwrite: bool = Field(False,
                            description="If True, updates the existing snapshot. If False, creates a new branch.")
