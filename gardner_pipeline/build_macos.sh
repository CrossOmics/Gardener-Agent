#!/bin/bash

# ==============================================================================
# gardener Pipeline macOS Build Script
#
# Usage:
#   1. Open Terminal
#   2. cd /path/to/gardener_pipeline
#   3. chmod +x build_macos.sh
#   4. ./build_macos.sh
# ==============================================================================

# Exit immediately if a command exits with a non-zero status
set -e

# Define Color Codes for output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}[INFO] Starting gardener macOS Build Process...${NC}"

# Define Root Directory (gardener_pipeline)
ROOT_DIR=$(pwd)

# Define Destination Directory for Python Binaries
RELEASE_DIR="$ROOT_DIR/gardener_frontend/electron/release"

# Ensure the release directory exists
mkdir -p "$RELEASE_DIR"

# ==============================================================================
# Step 0: Install Dependencies
# ==============================================================================
echo -e "${BLUE}[INFO] Installing Python dependencies...${NC}"

# Install requirements from the root directory
# Ensure you are in the correct virtual environment before running this script
pip install -r "$ROOT_DIR/requirements.txt"

echo -e "${GREEN}[SUCCESS] Dependencies installed.${NC}"


# ==============================================================================
# Step 1: Build Backend (gardenerBackend.app)
# ==============================================================================
echo -e "${BLUE}[INFO] Building Backend...${NC}"

# Navigate to backend directory
cd "$ROOT_DIR/gardener_backend"

# Clean previous PyInstaller builds to ensure fresh packaging
rm -rf build/ dist/

# Run PyInstaller using the macOS specific spec file
# --noconfirm overwrites output directory without asking
pyinstaller --clean --noconfirm run_server_macos.spec

# Check if the output exists
if [ -d "dist/gardenerBackend.app" ]; then
    echo -e "${GREEN}[SUCCESS] Backend built.${NC}"
else
    echo "Error: dist/gardenerBackend.app not found after build."
    exit 1
fi

# Clean destination if it exists to avoid move errors
rm -rf "$RELEASE_DIR/gardenerBackend.app"

# Move the generated .app to the electron release folder
mv "dist/gardenerBackend.app" "$RELEASE_DIR/"
echo -e "${GREEN}[SUCCESS] Moved gardenerBackend.app to electron/release/${NC}"


# ==============================================================================
# Step 2: Build Agent (gardenerAgent.app)
# ==============================================================================
echo -e "${BLUE}[INFO] Building Agent...${NC}"

# Navigate to agent directory
cd "$ROOT_DIR/gardener_agent"

# Clean previous PyInstaller builds
rm -rf build/ dist/

# Run PyInstaller using the macOS specific spec file
pyinstaller --clean --noconfirm run_agent_macos.spec

# Check if the output exists
if [ -d "dist/gardenerAgent.app" ]; then
    echo -e "${GREEN}[SUCCESS] Agent built.${NC}"
else
    echo "Error: dist/gardenerAgent.app not found after build."
    exit 1
fi

# Clean destination if it exists
rm -rf "$RELEASE_DIR/gardenerAgent.app"

# Move the generated .app to the electron release folder
mv "dist/gardenerAgent.app" "$RELEASE_DIR/"
echo -e "${GREEN}[SUCCESS] Moved gardenerAgent.app to electron/release/${NC}"


# ==============================================================================
# Step 3: Build Electron App (DMG)
# ==============================================================================
echo -e "${BLUE}[INFO] Building Electron Frontend...${NC}"

# Navigate to electron directory
cd "$ROOT_DIR/gardener_frontend/electron"

# Verify that the Python apps are in place before building
if [ ! -d "release/gardenerBackend.app" ] || [ ! -d "release/gardenerAgent.app" ]; then
    echo "Error: Python binaries are missing in release/ folder. Cannot proceed."
    exit 1
fi

# Detect CPU Architecture to run the correct npm command
ARCH=$(uname -m)

if [ "$ARCH" == "arm64" ]; then
    echo -e "${BLUE}[INFO] Detected Apple Silicon (M1/M2/M3). Running dist:mac:arm64...${NC}"
    npm run dist:mac:arm64
else
    echo -e "${BLUE}[INFO] Detected Intel Mac (x64). Running dist:mac...${NC}"
    npm run dist:mac
fi

# ==============================================================================
# Final Summary
# ==============================================================================
echo -e "${GREEN}====================================================${NC}"
echo -e "${GREEN} Build Complete! ${NC}"
echo -e "${GREEN}====================================================${NC}"
echo -e "Output DMG should be located in:"
echo -e "  $ROOT_DIR/gardener_frontend/electron/release/"