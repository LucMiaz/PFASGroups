"""
Pytest configuration file for PFASGroups test suite.

This file handles environment setup to avoid numpy floating-point exceptions
that can crash test collection on Windows.
"""
import os
import sys

# Set environment variables before importing numpy to prevent BLAS FPE crashes
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1' 
os.environ['NUMEXPR_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

# Disable floating-point exceptions at the OS level (Windows)
if sys.platform == 'win32':
    try:
        import ctypes
        # Disable floating-point exceptions in the Windows FPU control word
        # 0x9001F = enable all floating-point exceptions to be masked
        ctypes.windll.kernel32.SetErrorMode(0x0001 | 0x0002 | 0x8000)
        # Set FPU control word to mask all exceptions
        _control87 = ctypes.CDLL('msvcrt')._control87
        _control87.argtypes = [ctypes.c_uint, ctypes.c_uint]
        _control87.restype = ctypes.c_uint
        # MCW_EM = 0x0008001F (mask all exceptions)
        _control87(0x0008001F, 0x0008001F)
    except Exception:
        pass

# Try to disable numpy floating-point error checking early
try:
    import numpy as np
    # Disable floating-point exceptions
    np.seterr(all='ignore')
except ImportError:
    pass

import pytest
