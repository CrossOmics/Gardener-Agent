# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

block_cipher = None

# --- 1. Collect Metadata for Critical Libraries ---
# LangChain, LangGraph, and OpenAI integration often require package metadata
# to function correctly (e.g., for version checking or plugin discovery).
datas = []
try:
    datas += copy_metadata('langchain')
    datas += copy_metadata('langchain_core')
    datas += copy_metadata('langchain_community')
    datas += copy_metadata('langchain_openai')
    datas += copy_metadata('langgraph')
except Exception as e:
    print(f"Warning: Could not copy some metadata: {e}")

# --- 2. Define Hidden Imports ---
# Uvicorn loads these modules dynamically using import strings.
# PyInstaller cannot detect them automatically, so we list them here.
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
    # LangChain specific hidden imports
    'langchain_community.chat_models',
    'langchain_openai',
]

a = Analysis(
    ['run_agent.py'],
    pathex=['src'],              # CRITICAL: Adds 'src' to the import search path
    binaries=[],
    datas=datas,                 # Include the metadata we collected above
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['pytest', 'tests', 'tkinter'], # Exclude unnecessary testing libs
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

# --- 3. OneFile Configuration ---
# Bundle everything (scripts, binaries, data) into a single EXE file.
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,   # Include dynamic libraries (DLLs)
    a.zipfiles,   # Include compiled Python modules
    a.datas,      # Include non-code data (metadata, etc.)
    [],
    name='lotus-agent',          # The final output name (lotus-agent.exe)
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,                # Keep True for debugging errors
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    exclude_binaries=False,      # Must be False for OneFile mode
)