# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import copy_metadata

block_cipher = None

# --- 1. Collect Metadata ---
# LangChain and its plugins require metadata to be present at runtime.
datas = []
try:
    datas += copy_metadata('langchain')
    datas += copy_metadata('langchain_core')
    datas += copy_metadata('langchain_community')
    datas += copy_metadata('langchain_openai')
    datas += copy_metadata('langgraph')
except Exception as e:
    print(f"Warning: Could not copy metadata: {e}")

# --- 2. Hidden Imports ---
# Uvicorn loads these dynamically, so we must list them manually.
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
    'langchain_community.chat_models',
    'langchain_openai',
    'langchain',
    'langchain_core',
    'langchain_community',
    'langgraph',
    'planning', 
]

a = Analysis(
    ['run_agent.py'],
    pathex=['src'],  # Ensure 'src' is in the python path
    binaries=[],     # macOS usually finds dynamic libraries automatically
    datas=datas,
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

# --- 3. OneFile Configuration ---
# Bundle everything into a single executable file.
exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.zipfiles,
    a.datas,
    [],
    name='lotus-agent',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,  # Keep True for debugging
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    exclude_binaries=False, # Must be False for OneFile

)

# --- 4. macOS App Bundle ---
# Wraps the executable in a standard .app structure
app = BUNDLE(
    exe,
    name='LotusAgent.app',
    icon=None,
    bundle_identifier='com.bioagent.lotusagent',
    info_plist={
        'NSHighResolutionCapable': 'True',
        'LSBackgroundOnly': 'False',
    },
)