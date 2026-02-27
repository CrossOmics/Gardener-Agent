from typing import List, Union, Optional
from pydantic import Field
from dto.request.base_request import BaseAnalysisRequest


class MergeClustersRequest(BaseAnalysisRequest):
    """
    Request DTO for merging clusters.
    """
    method: str = Field("leiden", description="Clustering method column to modify, e.g., 'leiden', 'louvain'")

    clusters_to_merge: Union[List[str], List[int], List[float]] = Field(
        ...,
        description="List of cluster IDs to merge, e.g., ['1', '4', '8']"
    )

    new_label: Optional[str] = Field(
        None,
        description="Label for the merged group. If None, defaults to joined IDs (e.g., '1-4-8')."
    )

    overwrite: bool = Field(
        False,
        description="If True, updates the existing snapshot in-place. If False, creates a new snapshot branch."
    )
