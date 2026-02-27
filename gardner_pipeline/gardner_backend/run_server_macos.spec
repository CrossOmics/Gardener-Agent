# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata, collect_dynamic_libs
import sys
import os

block_cipher = None

# --- 1. Metadata Collection ---
metadata_packages = [
    'scikit-learn',
    'pandas',
    'numpy',
    'scipy',
    'scanpy',
    'anndata',
    'matplotlib',
    'seaborn',
    'cplearn',
]

all_metadata = []
for pkg in metadata_packages:
    try:
        all_metadata.extend(copy_metadata(pkg))
    except Exception as e:
        print(f"Warning: Could not copy metadata for {pkg}: {e}")

# --- 2. Hidden Imports ---
hidden_imports = [
    'uvicorn.logging',
    'uvicorn.loops',
    'uvicorn.loops.auto',
    'uvicorn.protocols',
    'uvicorn.protocols.http',
    'uvicorn.protocols.http.auto',
    'uvicorn.protocols.websockets',
    'uvicorn.protocols.websockets.auto',
    'uvicorn.lifespan.on',
    'engineio.async_drivers.asgi',
    'lotus',
    'lotus.core_analysis.cplearn'
]

# --- 3. Dynamic Binary Collection (The Fix) ---
# Finds TBB library dynamically for environment compatibility (CI/CD vs Local)
binaries = collect_dynamic_libs('tbb')

# --- 4. Custom Resource Inclusion ---
# Format: ( 'Source Path (relative to spec)', 'Dest Path (inside bundle)' )
# This ensures 'celltypist_models.json' and 'gseapy_libraries.json' are packaged.
project_resources = [
    (
        'src/infrastructure/database/resources',   # Source: Where your JSONs live
        'infrastructure/database/resources'        # Dest: Where Python code expects them
    )
]

a = Analysis(
    ['run_server.py'],
    pathex=['src'],
    # Use the dynamically collected binaries
    binaries=binaries,
    # ADDED: Merge metadata with your custom project resources
    datas=all_metadata + project_resources,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['pytest', 'tests', 'tkinter'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# --- 5. OneFile Configuration ---
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='lotus-backend',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    exclude_binaries=False,
)

# --- 6. macOS App Bundle ---
app = BUNDLE(
    exe,
    name='LotusBackend.app',
    icon=None,
    bundle_identifier='com.bioagent.lotusbackend',
    info_plist={
        'NSHighResolutionCapable': 'True',
        'LSBackgroundOnly': 'False',
    },
)