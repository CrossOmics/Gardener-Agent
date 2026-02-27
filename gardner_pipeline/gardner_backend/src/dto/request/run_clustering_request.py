from typing import Literal
from pydantic import Field

from dto.request.base_request import BaseAnalysisRequest


class RunClusteringRequest(BaseAnalysisRequest):

    # Clustering parameters
    method: Literal['leiden', 'louvain', 'cplearn'] = 'leiden'
    resolution: float = Field(0.5, description="List of resolutions to test")

    # Options
    run_hierarchical: bool = True
