<#
.SYNOPSIS
    gardener Pipeline Windows Build Script

.DESCRIPTION
    Automates the packaging of Python backends and the Electron frontend for Windows.

    Steps:
    0. Install Dependencies (requirements.txt + TBB)
    1. Clean and Build Backend (gardener-backend.exe)
    2. Clean and Build Agent (gardener-agent.exe)
    3. Move binaries to electron/release/
    4. Build Electron Installer (NSIS)

.NOTES
    Run this script from the 'gardener_pipeline' root directory.
    Ensure 'pyinstaller', 'npm', and 'pip' are in your system PATH.
#>

# Stop script on first error
$ErrorActionPreference = "Stop"

Write-Host "====================================================" -ForegroundColor Cyan
Write-Host "       Starting gardener Windows Build Process         " -ForegroundColor Cyan
Write-Host "====================================================" -ForegroundColor Cyan

# Define Root Directory (Current folder)
$RootDir = Get-Location

# Define Release Directory where binaries must go
$ReleaseDir = Join-Path $RootDir "gardener_frontend\electron\release"

# Create release directory if it doesn't exist
if (-not (Test-Path -Path $ReleaseDir)) {
    New-Item -ItemType Directory -Path $ReleaseDir | Out-Null
    Write-Host "[INFO] Created release directory." -ForegroundColor Gray
}

# ==============================================================================
# Step 0: Install Dependencies
# ==============================================================================
Write-Host "`n[STEP 0] Installing Python dependencies..." -ForegroundColor Yellow

# Check if requirements.txt exists
$ReqFile = Join-Path $RootDir "requirements.txt"

if (Test-Path $ReqFile) {
    # Install dependencies from file
    Write-Host "  -> Installing from requirements.txt..." -ForegroundColor Gray
    pip install -r $ReqFile

    # Install TBB specifically for Windows Numba support
    Write-Host "  -> Installing TBB (Required for Numba on Windows)..." -ForegroundColor Gray
    pip install tbb

    if ($LASTEXITCODE -eq 0) {
        Write-Host "[SUCCESS] Dependencies installed." -ForegroundColor Green
    } else {
        Write-Error "Failed to install dependencies. Please check your network or python environment."
    }
} else {
    Write-Warning "requirements.txt not found in root directory. Skipping installation."
}

# ==============================================================================
# Step 1: Build Backend (gardener-backend.exe)
# ==============================================================================
Write-Host "`n[STEP 1] Building Backend..." -ForegroundColor Yellow

# Navigate to backend directory
Set-Location -Path (Join-Path $RootDir "gardener_backend")

# Clean previous PyInstaller builds
if (Test-Path "build") { Remove-Item -Path "build" -Recurse -Force }
if (Test-Path "dist") { Remove-Item -Path "dist" -Recurse -Force }

# Run PyInstaller
# --clean: Clean PyInstaller cache
# --noconfirm: Overwrite output directory without asking
Write-Host "  -> Running PyInstaller for run_server.spec..." -ForegroundColor Gray
pyinstaller --clean --noconfirm run_server.spec
# py -3.11 -m PyInstaller --clean --noconfirm run_server.spec
# Verify output exists
$BackendExe = "dist\gardener-backend.exe"
if (-not (Test-Path $BackendExe)) {
    Write-Error "Build Failed: gardener-backend.exe not found."
}

# Move to release folder (overwrite if exists)
Move-Item -Path $BackendExe -Destination $ReleaseDir -Force
Write-Host "[SUCCESS] Moved gardener-backend.exe to electron/release/" -ForegroundColor Green

# ==============================================================================
# Step 2: Build Agent (gardener-agent.exe)
# ==============================================================================
Write-Host "`n[STEP 2] Building Agent..." -ForegroundColor Yellow

# Navigate to agent directory
Set-Location -Path (Join-Path $RootDir "gardener_agent")

# Clean previous PyInstaller builds
if (Test-Path "build") { Remove-Item -Path "build" -Recurse -Force }
if (Test-Path "dist") { Remove-Item -Path "dist" -Recurse -Force }

# Run PyInstaller
Write-Host "  -> Running PyInstaller for run_agent.spec..." -ForegroundColor Gray
pyinstaller --clean --noconfirm run_agent.spec

# Verify output exists
$AgentExe = "dist\gardener-agent.exe"
if (-not (Test-Path $AgentExe)) {
    Write-Error "Build Failed: gardener-agent.exe not found."
}

# Move to release folder (overwrite if exists)
Move-Item -Path $AgentExe -Destination $ReleaseDir -Force
Write-Host "[SUCCESS] Moved gardener-agent.exe to electron/release/" -ForegroundColor Green

# ==============================================================================
# Step 3: Build Electron App (NSIS Installer)
# ==============================================================================
Write-Host "`n[STEP 3] Building Electron Frontend (NSIS)..." -ForegroundColor Yellow

# Navigate to electron directory
Set-Location -Path (Join-Path $RootDir "gardener_frontend\electron")

# Verify binaries exist before starting electron-builder
if (-not (Test-Path "release\gardener-backend.exe") -or -not (Test-Path "release\gardener-agent.exe")) {
    Write-Error "Missing Python binaries in release/ folder. Cannot proceed."
}

# Install npm dependencies (Optional, typically needed once)
# Write-Host "  -> Installing npm dependencies..."
# npm install

# Build Windows Installer
Write-Host "  -> Running npm run dist:win..." -ForegroundColor Gray
# Using cmd /c to ensure npm runs correctly in PowerShell execution environments
cmd /c "npm run dist:win"

if ($LASTEXITCODE -eq 0) {
    Write-Host "`n====================================================" -ForegroundColor Green
    Write-Host "             Windows Build Complete!                " -ForegroundColor Green
    Write-Host "====================================================" -ForegroundColor Green
    Write-Host "Installer location:"
    Write-Host "$ReleaseDir" -ForegroundColor White
} else {
    Write-Error "Electron build failed."
}

# Return to root
Set-Location $RootDir