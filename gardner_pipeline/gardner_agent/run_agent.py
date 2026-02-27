import multiprocessing
import sys
import os
import uvicorn

# 1. Standard boilerplate for Windows to prevent infinite subprocess spawning
if __name__ == "__main__":
    multiprocessing.freeze_support()

    # 2. Add 'src' to the system path.
    #    This allows us to do "from main import app" even though main.py is inside src/.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    src_path = os.path.join(current_dir, "src")
    sys.path.insert(0, src_path)

    # 3. Import the FastAPI app.
    #    We import it inside the main block to avoid recursive import issues.
    try:
        from main import app
    except ImportError as e:
        # Debugging block: helps you see why the import failed in the compiled EXE
        print(f"CRITICAL ERROR: Could not import main app: {e}")
        import traceback
        traceback.print_exc()
        input("Press Enter to exit...")
        sys.exit(1)

    print("Starting Agent Server on port 41889...")

    # 4. Start Uvicorn
    #    IMPORTANT: Pass the 'app' object directly, not the string "main:app".
    #    'reload' must be False for packaged apps.
    uvicorn.run(
        app,
        host="0.0.0.0",
        port=41889,
        reload=False,
        workers=1,
        log_level="info"
    )