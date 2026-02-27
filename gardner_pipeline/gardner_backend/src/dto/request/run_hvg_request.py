from pydantic import Field
from typing import Literal

from dto.request.base_request import BaseAnalysisRequest


class RunHVGRequest(BaseAnalysisRequest):

    n_top_genes: int = Field(2000, ge=500, le=10000, description="Target number of HVGs")
    flavor: Literal['seurat', 'cell_ranger', 'seurat_v3'] = 'seurat'
    target_sum: float = 1e4
