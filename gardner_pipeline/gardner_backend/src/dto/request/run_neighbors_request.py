from dto.request.base_request import BaseAnalysisRequest


class RunNeighborsRequest(BaseAnalysisRequest):
    n_neighbors: int = 15
    n_pcs: int = 30
