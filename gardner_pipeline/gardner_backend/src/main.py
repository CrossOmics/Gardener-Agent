
import uvicorn
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from contextlib import asynccontextmanager

from core.task_executor import task_executor
from infrastructure.database.connection import get_default_db_manager
from infrastructure.database.database_setup import db_setup
from infrastructure.filesystem.logging import setup_logging
from infrastructure.workspace_context import workspace_path_manager

from controller.agent_chat_controller import router as agent_chat_router
from controller.annotation_controller import router as annotation_router
from controller.clustering_controller import router as clustering_router
from controller.dge_controller import router as dge_router
from controller.dim_reduction_controller import router as dim_reduction_router
from controller.preprocessing_controller import router as preprocessing_router
from controller.project_management_controller import router as project_management_router
from controller.run_whole_pipeline_controller import router as pipeline_router
from controller.search_annotation_methods_controller import router as search_methods_router
from controller.snapshot_controller import router as snapshot_router
from controller.user_preference_controller import router as user_preference_router


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application startup and shutdown lifecycle."""

    print("System starting up...")
    # Initialize workspace context
    workspace_path_manager.initialize()
    # Initialize logging
    setup_logging(debug_mode=True)

    # Initialize database connections
    db_setup(ensure_schema=True)

    yield

    print("System shutting down...")
    # Shutdown thread pool
    task_executor.shutdown()
    # Close database connection
    get_default_db_manager().close_connection()


app = FastAPI(
    title="Gardener Biological Agent API",
    description="Backend engine for single-cell analysis",
    version="0.1.0",
    lifespan=lifespan,
)

# Allow all origins for local desktop / WebView usage
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Manually include all controllers to expose API endpoints.

# 1. Core Project & Data Management
app.include_router(project_management_router)
app.include_router(user_preference_router)
app.include_router(snapshot_router)

# 2. Analysis Pipeline Steps
app.include_router(preprocessing_router)
app.include_router(dim_reduction_router)
app.include_router(clustering_router)
app.include_router(dge_router)
app.include_router(annotation_router)

# 3. Automation & Utilities
app.include_router(pipeline_router)
app.include_router(search_methods_router)

# 4. AI Agent
app.include_router(agent_chat_router)


@app.get("/health")
def health_check():
    """Health check endpoint for frontend readiness detection."""
    return {"status": "ok", "version": "0.1.0"}


if __name__ == "__main__":
    # For local development only.
    # In production, the server is started by PyWebView or an external launcher.
    uvicorn.run(
        "main:app",
        host="127.0.0.1",
        port=41888,
        reload=True,
        access_log=False
    )
