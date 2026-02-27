# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import (
    copy_metadata,
    collect_dynamic_libs,
    collect_submodules,
    collect_data_files
)
import sys
import os
import glob

block_cipher = None

# --- 1. Collect Package Metadata ---
metadata_packages = [
    'scikit-learn',
    'scanpy',
    'anndata',
    'pandas',
    'numpy',
    'scipy',
    'matplotlib',
    'seaborn',
    'statsmodels',
    'networkx',
    'tbb',
    'cplearn',
]

all_datas = []
for pkg in metadata_packages:
    try:
        all_datas.extend(copy_metadata(pkg))
    except Exception as e:
        print(f"Warning: Could not copy metadata for {pkg}: {e}")

# Collect data files for sklearn and scipy
try:
    all_datas.extend(collect_data_files('sklearn'))
    all_datas.extend(collect_data_files('scipy'))
except Exception as e:
    print(f"Warning: Could not collect data files: {e}")


# --- 2. Hidden Imports ---
hidden_imports = [
    # Uvicorn async server
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

    # Scikit-learn Cython extensions
    'sklearn.neighbors._partition_nodes',
    'sklearn.utils._cython_blas',
    'sklearn.utils._typedefs',
    'sklearn.utils._weight_vector',
    'sklearn.tree._utils',
    'sklearn.neighbors._quad_tree',
    'sklearn.ensemble._gradient_boosting',

    # SciPy Cython extensions
    'scipy.special.cython_special',
    'scipy.linalg.cython_blas',
    'scipy.linalg.cython_lapack',
    'scipy.sparse.csgraph._validation',

    # Other dependencies
    'pkg_resources.py2_warn',
    'lotus',
    'lotus.core_analysis.cplearn',
    'cplearn',
]

# Auto-collect sklearn submodules
try:
    hidden_imports.extend(collect_submodules('sklearn.utils'))
    hidden_imports.extend(collect_submodules('sklearn.metrics'))
except Exception as e:
    print(f"Warning: Could not collect sklearn submodules: {e}")


# --- 3. Project Resources ---
project_resources = [
    (
        'src/infrastructure/database/resources',
        'infrastructure/database/resources'
    )
]


# --- 4. Collect Dynamic Libraries (DLLs) ---
binaries = []

# Auto-collect binaries from key packages
binary_packages = ['tbb', 'numpy', 'scipy', 'sklearn', 'cplearn']
for pkg in binary_packages:
    try:
        print(f"Collecting binaries for {pkg}...")
        collected = collect_dynamic_libs(pkg)
        binaries.extend(collected)
        print(f"  Collected {len(collected)} files")
    except Exception as e:
        print(f"  Warning: {e}")

# Windows-specific: Scan Conda Library/bin for TBB and MKL DLLs
if sys.platform == 'win32':
    conda_bin_path = os.path.join(sys.prefix, 'Library', 'bin')

    if os.path.exists(conda_bin_path):
        print(f"Scanning Conda bin: {conda_bin_path}")

        # Collect TBB and MKL DLLs (numpy dependency)
        dll_patterns = ['tbb*.dll', 'mkl*.dll']

        for pattern in dll_patterns:
            for dll_path in glob.glob(os.path.join(conda_bin_path, pattern)):
                dll_name = os.path.basename(dll_path)
                # Avoid duplicates
                if not any(dll_name in b[0] for b in binaries):
                    print(f"  Found: {dll_name}")
                    binaries.append((dll_path, '.'))

    # Fallback: Explicitly add critical DLLs if not found
    critical_dlls = ['tbb12.dll', 'tbbmalloc.dll', 'tbbmalloc_proxy.dll']

    for dll_name in critical_dlls:
        if not any(dll_name in b[0] for b in binaries):
            possible_paths = [
                os.path.join(sys.prefix, 'Library', 'bin', dll_name),
                os.path.join(sys.prefix, 'DLLs', dll_name),
                os.path.join(sys.prefix, dll_name),
            ]

            for path in possible_paths:
                if os.path.exists(path):
                    print(f"  Fallback: {dll_name}")
                    binaries.append((path, '.'))
                    break

print(f"\nTotal binaries: {len(binaries)}")


# --- 5. Analysis Phase ---
a = Analysis(
    ['run_server.py'],
    pathex=['src', '.'],
    binaries=binaries,
    datas=all_datas + project_resources,
    hiddenimports=hidden_imports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['pytest', 'tests', 'tkinter', 'matplotlib.tests'],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)


# --- 6. OneFile Executable ---
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
    upx=False,  # Disable UPX compression (critical for DLL loading in OneFile)
    upx_exclude=[],
    runtime_tmpdir=None,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    exclude_binaries=False,  # OneFile mode: embed all binaries
)