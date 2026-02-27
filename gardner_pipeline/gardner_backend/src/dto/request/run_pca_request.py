from dto.request.base_request import BaseAnalysisRequest


class RunPCARequest(BaseAnalysisRequest):
    n_comps: int = 50
    svd_solver: str = 'arpack'