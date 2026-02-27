import multiprocessing
import sys
import os
import platform

# Platform Detection
IS_MACOS = platform.system() == 'Darwin'
IS_WINDOWS = platform.system() == 'Windows'

# Force joblib to use threading instead of spawning processes
os.environ["LOKY_MAX_CPU_COUNT"] = "4"  # Limit threads
os.environ["LOKY_PICKLER"] = "pickle"   # Use standard pickle

# Limit thread counts to prevent resource exhaustion
os.environ["OMP_NUM_THREADS"] = "4"
os.environ["MKL_NUM_THREADS"] = "4"
os.environ["OPENBLAS_NUM_THREADS"] = "4"

# CRITICAL: Tell joblib to prefer threading over multiprocessing
# This must be set BEFORE importing any joblib-using libraries
import joblib

joblib.parallel.BACKENDS['loky'] = joblib.parallel.BACKENDS['threading']

# === macOS-specific fixes ===
if IS_MACOS:
    # Fix for macOS multiprocessing in frozen app
    multiprocessing.set_start_method('spawn', force=True)

if __name__ == "__main__":
    # Windows/macOS PyInstaller requirement
    if IS_WINDOWS or IS_MACOS:
        multiprocessing.freeze_support()

    current_dir = os.path.dirname(os.path.abspath(__file__))
    src_path = os.path.join(current_dir, "src")
    sys.path.insert(0, src_path)

    try:
        import uvicorn
        from src.main import app
    except ImportError as e:
        print(f"Critical Import Error: {e}")
        import traceback

        traceback.print_exc()

        if not IS_MACOS:
            input("Press Enter to exit...")
        sys.exit(1)

    print("Starting Backend Service on port 41888...")

    uvicorn.run(
        app,
        host="127.0.0.1",
        port=41888,
        reload=False,
        workers=1,
        log_level="info"
    )