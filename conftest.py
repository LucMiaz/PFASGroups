"""
Pytest configuration file for HalogenGroups test suite.

Note: There is a known issue with Python 3.14 + numpy on Windows where
numpy's BLAS FPE check causes crashes during test collection. 

Workaround: Run individual test files directly with Python instead of pytest:
    python -m tests.test_pfas_identification
"""
import pytest
