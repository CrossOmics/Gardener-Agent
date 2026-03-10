import uvicorn
import sys
import io
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware

from controller.agent_chat_controller import router as agent_chat_router
from controller.agent_annotation_controller import router as annotation_router

# Force UTF-8 Encoding for Windows Console
if sys.stdout:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
if sys.stderr:
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

app = FastAPI(title="AI Agent Server")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Register Routers
app.include_router(agent_chat_router)
app.include_router(annotation_router)

if __name__ == "__main__":
    # Start the Uvicorn server on port 41889
    # reload=True is good for dev, consider removing for production
    uvicorn.run("main:app", host="0.0.0.0", port=41889, reload=True)
