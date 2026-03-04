"""
Pytest configuration file for HalogenGroups test suite.

Note: There is a known issue with Python 3.14 + numpy on Windows where
numpy's BLAS FPE check causes crashes during test collection because the
conda env's DLL directories are not added to the Windows DLL search path
when running pytest directly (outside of `conda run`).

The fix below patches os.add_dll_directory() and PATH before any imports.
"""
import os
import sys

# Fix for numpy/rdkit DLL load failures on Python 3.14 + Windows (conda envs).
# When pytest is launched without `conda run`, the conda env's Library/bin is
# not in PATH, so BLAS and other DLLs cannot be found.
_conda_prefix = os.path.dirname(os.path.dirname(sys.executable))
_dll_dirs = [
    os.path.join(_conda_prefix, "Library", "bin"),
    os.path.join(_conda_prefix, "Library", "mingw-w64", "bin"),
    os.path.join(_conda_prefix, "Library", "usr", "bin"),
    os.path.join(_conda_prefix, "bin"),
]
for _d in _dll_dirs:
    if os.path.isdir(_d):
        if hasattr(os, "add_dll_directory"):
            try:
                os.add_dll_directory(_d)
            except OSError:
                pass
        if _d not in os.environ.get("PATH", ""):
            os.environ["PATH"] = _d + os.pathsep + os.environ.get("PATH", "")

# Suppress duplicate OpenMP library warnings (common with conda rdkit + numpy)
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import pytest
