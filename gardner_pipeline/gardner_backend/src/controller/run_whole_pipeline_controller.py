
from typing import Dict, List
from fastapi import APIRouter, Depends, HTTPException, status
from loguru import logger

from core.task_executor import task_executor
from dto.response.run_whole_pipeline_response import RunWholePipelineResponse

from service.preprocessing_service import PreprocessingService
from service.clustering_service import ClusteringService
from service.dge_service import DGEService
from service.annotation_service import AnnotationService

from dto.request.run_whole_pipeline_request import RunWholePipelineRequest

from util.id_generate_utils import generate_business_id

router = APIRouter(
    prefix="/api/v1/pipeline",
    tags=["Whole Pipeline Automation"]
)


def _run_pipeline_task(
        request: RunWholePipelineRequest,
        prep_service: PreprocessingService,
        cluster_service: ClusteringService,
        deg_service: DGEService,
        annot_service: AnnotationService
) -> RunWholePipelineResponse:
    """Synchronous function to be executed in a thread."""
    pipeline_run_id = generate_business_id("run_pipe")
    snapshots_map: Dict[str, str] = {}
    steps_completed: List[str] = []

    pid = request.project_id
    did = request.dataset_id

    try:
        # Step 1: Preprocessing
        logger.info(f"[{pipeline_run_id}] State 1: Full Preprocessing")
        pp_req = request.preprocessing_params
        pp_req.project_id = pid
        pp_req.dataset_id = did
        pp_res = prep_service.full_preprocessing(pp_req)
        snapshots_map["preprocessing"] = getattr(pp_res, 'snapshot_id', 'done')
        steps_completed.extend(pp_res.executed_steps)

        # Step 2: Clustering
        logger.info(f"[{pipeline_run_id}] State 2: Clustering")
        cluster_req = request.clustering_params
        cluster_req.project_id = pid
        cluster_req.dataset_id = did
        cluster_res = cluster_service.run_clustering(cluster_req)
        snapshots_map["clustering"] = cluster_res.snapshot_id
        steps_completed.append("Clustering")

        # Step 3: DGE
        logger.info(f"[{pipeline_run_id}] State 3: DEG")
        deg_req = request.deg_params
        deg_req.project_id = pid
        deg_req.dataset_id = did
        deg_res = deg_service.run_dge_analysis(deg_req)
        snapshots_map["deg"] = deg_res.snapshot_id
        steps_completed.append("DEG")

        # Step 4: Annotation
        logger.info(f"[{pipeline_run_id}] State 4: Full Annotation")
        annot_req = request.annotation_params
        annot_req.project_id = pid
        annot_req.dataset_id = did
        full_annot_res = annot_service.run_full_annotation(annot_req)
        snapshots_map["annotation_full"] = full_annot_res.snapshot_id
        steps_completed.append("Full Annotation")

        logger.success(f"Pipeline {pipeline_run_id} Completed.")
        return RunWholePipelineResponse(
            pipeline_run_id=pipeline_run_id,
            final_snapshot_id=full_annot_res.snapshot_id,
            status="SUCCESS",
            msg="Whole pipeline executed successfully.",
            steps_completed=steps_completed,
            snapshots=snapshots_map
        )
    except Exception as e:
        logger.error(f"Pipeline {pipeline_run_id} Failed: {e}")
        # Re-raise as HTTPException to be handled by FastAPI
        raise HTTPException(status_code=500, detail=str(e))


@router.post("/run", response_model=RunWholePipelineResponse, status_code=status.HTTP_201_CREATED)
async def run_whole_pipeline(
        request: RunWholePipelineRequest,
        prep_service: PreprocessingService = Depends(),
        cluster_service: ClusteringService = Depends(),
        deg_service: DGEService = Depends(),
        annot_service: AnnotationService = Depends()
):
    """
    Asynchronously trigger the whole pipeline execution in a separate thread.
    """
    return await task_executor.run_in_thread(
        _run_pipeline_task,
        request,
        prep_service,
        cluster_service,
        deg_service,
        annot_service
    )
