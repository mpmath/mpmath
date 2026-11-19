import sys

import pytest

import mpmath


def pytest_report_header(config):
    print(f"mpmath backend: {mpmath.libmp.backend.BACKEND}")
    print(f"mpmath mp class: {mpmath.mp!r}")
    print(f"mpmath version: {mpmath.__version__}")
    print(f"Python version: {sys.version}")


def pytest_configure(config):
    config.addinivalue_line('markers', 'slow: marks tests as slow')


@pytest.fixture(autouse=True)
def reset_mp_globals():
    mpmath.mp.prec = sys.float_info.mant_dig
    mpmath.mp.pretty = False
    mpmath.mp.rounding = 'n'
    mpmath.mp.pretty_dps = "str"
    mpmath.mp.shortest_str = False
    mpmath.iv.prec = mpmath.mp.prec
    mpmath.iv.pretty = False
