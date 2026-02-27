from typing import Annotated
from fastapi import APIRouter, Depends, HTTPException, status, Query
from loguru import logger

from core.task_executor import task_executor
from dto.request.display_name_request import DisplayNameRequest
# Import DTOs
from dto.request.project_management_request import CreateDatasetRequest, CreateProjectRequest
from dto.request.rename_request import RenameRequest
from dto.response.display_name_response import DisplayNameResponse
from dto.response.project_response import ProjectResponse, ProjectEntityMapResponse
from dto.response.rename_response import RenameResponse
from infrastructure.workspace_context import workspace_path_manager
from service.agent_chat_service import AgentChatService

# Import Services
from service.dataset_service import DatasetService
from service.project_management_service import ProjectManagementService
from service.snapshot_service import SnapshotService

router = APIRouter(prefix="/api/v1/project", tags=["Project Management"])

# Dependency Injection Setup
ProjectServiceDep = Annotated[ProjectManagementService, Depends()]
DatasetServiceDep = Annotated[DatasetService, Depends()]
ChatServiceDep = Annotated[AgentChatService, Depends()]
SnapshotServiceDep = Annotated[SnapshotService, Depends()]


@router.post(
    "/create",
    response_model=ProjectResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Create New Project",
    description="Creates a new empty project workspace with metadata."
)
async def create_new_project(
        request: CreateProjectRequest,
        project_service: ProjectServiceDep,
        chat_service: ChatServiceDep
):
    """
    Creates a new empty project structure and initializes a default chat session.
    """
    def _task():
        logger.info(f"Creating new project: {request.project_name}")

        # 1. Create Project Metadata & Workspace
        new_project = project_service.create_project(
            project_name=request.project_name,
            description=request.description
        )

        if not new_project:
            raise HTTPException(
                status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                detail="Failed to create project metadata."
            )

        logger.success(f"Project created successfully: {new_project.project_id}")

        # 2. Initialize Default Chat Session
        # Default to empty string if session creation fails
        sid = ""
        try:
            default_session = chat_service.start_new_session(
                project_id=new_project.project_id,
                initial_name="Default Session"
            )
            sid = default_session.get("session_id", "")
            logger.info(f"Initialized default session: {sid}")
        except Exception as e:
            # We log the error but allow project creation to succeed
            logger.warning(f"Failed to initialize default session for {new_project.project_id}: {e}")

        # 3. Construct Response
        return ProjectResponse(
            project_id=new_project.project_id,
            project_name=new_project.project_name,
            workspace_path=new_project.project_path,
            session_id=sid
        )

    try:
        return await task_executor.run_in_thread(_task)
    except Exception as e:
        logger.exception("Unexpected error during project creation")
        raise HTTPException(status_code=500, detail=str(e))


@router.post(
    "/dataset/import",
    response_model=ProjectResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Import Dataset to Project",
    description="Ingests a raw dataset file from a local path into an existing project."
)
async def create_new_dataset(
        request: CreateDatasetRequest,
        project_service: ProjectServiceDep,
        dataset_service: DatasetServiceDep
):
    """
    Adds a new dataset to an existing project.

    Workflow:
    1. Verify Project Exists (using provided project_id).
    2. Validate Local File.
    3. Convert & Save standardized .h5ad file.
    4. Create Dataset Record with biological metadata.
    """
    def _task():
        logger.info(f"Importing dataset to Project {request.project_id} from {request.local_file_path}")

        # Retrieve Existing Project
        # need the project object to ensure it exists and to get the path
        target_project = project_service.get_project_by_id(request.project_id)

        if not target_project:
            logger.warning(f"Project ID {request.project_id} not found.")
            raise HTTPException(
                status_code=status.HTTP_404_NOT_FOUND,
                detail=f"Project {request.project_id} does not exist."
            )

        # Import Dataset (IO + DB)
        # This service method handles file copying, H5AD conversion, and DB insertion
        new_dataset = dataset_service.import_dataset_from_local(
            local_file_path=request.local_file_path,
            project=target_project,
            dataset_name=request.dataset_name,
            organism=request.organism,  # Biological meta passed here
            tissue_type=request.tissue_type,
            details=request.description  # Map description to 'details' field
        )

        if not new_dataset:
            raise HTTPException(
                status_code=status.HTTP_400_BAD_REQUEST,
                detail="Failed to import dataset. Check file integrity or format."
            )

        logger.success(f"Dataset {new_dataset.dataset_id} imported to {target_project.project_name}")

        # Return Combined Response
        return ProjectResponse(
            project_id=target_project.project_id,
            project_name=target_project.project_name,
            dataset_id=new_dataset.dataset_id,
            dataset_name=new_dataset.dataset_name,
            workspace_path=target_project.project_path
        )

    try:
        return await task_executor.run_in_thread(_task)
    except FileNotFoundError as e:
        logger.error(f"File error: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        logger.exception("System error during dataset import")
        raise HTTPException(status_code=500, detail=str(e))


@router.get(
    "/root",
    response_model=str,
    status_code=status.HTTP_200_OK,
    summary="Get Workspace Root Path",
    description="Returns the absolute file system path of the current workspace root."
)
async def get_workspace_root_path():
    """
    Retrieves the configured absolute path for the data root.
    Useful for frontend or users to know where files are stored locally.
    """
    try:
        root_path = str(workspace_path_manager.root.resolve())
        logger.info(f"Requested workspace root: {root_path}")
        return root_path
    except Exception as e:
        logger.error(f"Failed to retrieve workspace root: {e}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="Could not determine workspace root."
        )


@router.delete(
    "/manage/delete/{project_id}",
    status_code=status.HTTP_200_OK,
    summary="Delete Project",
    description="Deletes a project record and its associated workspace files."
)
async def delete_project(
        project_id: str,
        project_service: ProjectServiceDep
):
    """
    Deletes a project.

    1. Validates project existence.
    2. Removes physical files from storage.
    3. Deletes the database record.
    """
    def _task():
        logger.info(f"Deleting project: {project_id}")
        return project_service.delete_project_by_id(project_id)

    try:
        return await task_executor.run_in_thread(_task)
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        logger.exception(f"Failed to delete project {project_id}")
        raise HTTPException(status_code=500, detail=str(e))


@router.delete(
    "/dataset/delete",
    status_code=status.HTTP_200_OK,
    summary="Delete Dataset",
    description="Deletes a dataset record and its associated files. Optionally retains the latest snapshot."
)
async def delete_dataset(
        project_service: ProjectServiceDep,
        dataset_id: str = Query(..., description="The unique business ID of the dataset to delete"),
        keep_final: bool = Query(False,
                                 description="If true, the latest annotation snapshot of this dataset will be preserved"),
):
    """
    Deletes a dataset based on query parameters.

    1. Validates dataset existence via dataset_id.
    2. Removes physical files (raw data & snapshots).
    3. If keep_latest is True, skips the most recent snapshot during cleanup.
    4. Deletes the database record.
    """
    def _task():
        logger.info(f"Deleting dataset: {dataset_id} (keep_latest={keep_final})")
        return project_service.delete_dataset(dataset_id, keep_final=keep_final)

    try:
        return await task_executor.run_in_thread(_task)
    except ValueError as e:
        logger.error(f"Validation error: {e}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        logger.exception(f"Failed to delete dataset {dataset_id}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"An error occurred while deleting the dataset: {str(e)}"
        )


@router.patch(
    "/rename",
    response_model=RenameResponse,
    status_code=status.HTTP_200_OK,
    summary="Update Display Name",
    description="Updates the user-friendly display name of a Project, Dataset, or Snapshot based on its system ID."
)
async def update_display_name(
        request: RenameRequest,
        project_service: ProjectServiceDep,
        dataset_service: DatasetServiceDep,
        snapshot_service: SnapshotServiceDep
):
    """
    Updates the display name of a specific entity.
    Routes the request to the appropriate service based on 'id_type'.
    """
    def _task():
        logger.info(f"Renaming {request.id_type} '{request.current_id}' to '{request.new_name}'")

        if request.id_type == 'snapshot':
            # Calls the method defined in SnapshotService (handles its own verification)
            result = snapshot_service.update_snapshot_name(
                snapshot_id=request.current_id,
                new_name=request.new_name
            )

        elif request.id_type == 'dataset':
            # Assuming DatasetService has a similar update_dataset_name method
            result = dataset_service.update_dataset_name(
                dataset_id=request.current_id,
                new_name=request.new_name
            )

        elif request.id_type == 'project':
            # Assuming ProjectManagementService has an update_project_name method
            result = project_service.update_project_name(
                project_id=request.current_id,
                new_name=request.new_name
            )

        else:
            # Pydantic Literal validation should prevent this, but added for safety
            raise ValueError(f"Invalid id_type: {request.id_type}")

        # Construct unified response
        return RenameResponse(
            msg=f"{request.id_type.capitalize()} renamed successfully.",
            id_type=request.id_type,
            current_id=request.current_id,
            new_name=request.new_name
        )

    try:
        return await task_executor.run_in_thread(_task)
    except ValueError as e:
        logger.warning(f"Validation error during rename: {e}")
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        logger.exception(f"System error renaming {request.id_type} {request.current_id}")
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
                            detail="An error occurred while renaming the entity.")


@router.post(
    "/display/names",
    status_code=status.HTTP_200_OK,
    summary="Get Display Name",
    description="Retrieves the user-friendly display name of a Project, Dataset, or Snapshot based on its system ID."
)
async def get_display_name(
        request: DisplayNameRequest,
        project_service: ProjectServiceDep,
        dataset_service: DatasetServiceDep,
        snapshot_service: SnapshotServiceDep
):
    """
    Retrieves the display name of a specific entity based on its ID and Type.
    """
    def _task():
        id_type = request.id_type
        ids = request.current_ids
        logger.info(f"Bulk fetching names for {len(ids)} {id_type}(s)")

        results_map = {}

        # Route to the appropriate service based on type
        if id_type == 'snapshot':
            # Resolve names for each snapshot ID
            for sid in ids:
                try:
                    results_map[sid] = snapshot_service.get_snapshot_name(sid)
                except Exception:
                    results_map[sid] = ""

        elif id_type == 'dataset':
            for did in ids:
                try:
                    results_map[did] = dataset_service.get_dataset_name(did)
                except Exception:
                    results_map[did] = ""

        elif id_type == 'project':
            for pid in ids:
                try:
                    results_map[pid] = project_service.get_project_name(pid)
                except Exception:
                    results_map[pid] = ""

        return DisplayNameResponse(
            id_type=id_type,
            names=results_map
        )

    try:
        return await task_executor.run_in_thread(_task)
    except Exception as e:
        logger.exception(f"System error during bulk name retrieval for {request.id_type}")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail="An error occurred while retrieving the entity names."
        )


@router.get(
    "/display/map/{project_id}",
    response_model=ProjectEntityMapResponse,
    status_code=status.HTTP_200_OK,
    summary="Get Project Entity Name Map",
    description="Retrieves a full hierarchical map of IDs and Display Names (Project -> Datasets -> Snapshots) for a specific project."
)
async def get_project_entity_map(
        project_id: str,
        project_service: ProjectServiceDep,
):
    """
    Endpoint to fetch the nested ID-to-Name mapping for an entire project.
    """
    def _task():
        logger.info(f"Retrieving global entity map for project '{project_id}'")
        return project_service.get_project_entity_map(project_id)

    try:
        return await task_executor.run_in_thread(_task)
    except ValueError as e:
        logger.warning(f"Project not found: {e}")
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail=str(e))
    except Exception as e:
        logger.exception(f"System error generating entity map for project '{project_id}'")
        raise HTTPException(
            status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
            detail=f"Failed to retrieve entity map for project '{project_id}'"
        )
