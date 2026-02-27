from typing import Literal, Optional
from pydantic import BaseModel, Field

from dto.request.base_request import BaseAnalysisRequest


class RunDGERequest(BaseAnalysisRequest):
    # Analysis Parameters
    groupby: str = Field('leiden', description="Column in adata.obs to group by (usually clustering result)")
    method: Literal['wilcoxon', 't-test', 'logreg'] = Field('wilcoxon', description="Statistical test method")
    n_top_genes: int = Field(25, ge=5, le=100, description="Number of top genes to visualize")
    use_raw: bool = Field(True, description="Whether to use adata.raw for calculation")
