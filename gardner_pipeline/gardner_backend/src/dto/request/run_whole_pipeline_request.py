from typing import Optional
from pydantic import Field
from dto.request.base_request import BaseAnalysisRequest

from dto.request.full_preprocessing_request import FullPreprocessingRequest
from dto.request.run_clustering_request import RunClusteringRequest
from dto.request.run_dge_request import RunDGERequest
from dto.request.annotation_request import RunFullAnnotationRequest


class RunWholePipelineRequest(BaseAnalysisRequest):
    # 1. Preprocessing
    preprocessing_params: Optional[FullPreprocessingRequest] = Field(
        default_factory=FullPreprocessingRequest
    )

    # 2. Clustering
    clustering_params: RunClusteringRequest

    # 3. DEG
    deg_params: Optional[RunDGERequest] = Field(
        default_factory=RunDGERequest
    )

    # 4. Annotation
    annotation_params: Optional[RunFullAnnotationRequest] = Field(
        default_factory=RunFullAnnotationRequest
    )
