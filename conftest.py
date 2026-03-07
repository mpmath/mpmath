import sys

import pytest

import mpmath


def pytest_report_header(config):
    print("mpmath backend: %s" % mpmath.libmp.backend.BACKEND)
    print("mpmath mp class: %s" % repr(mpmath.mp))
    print("mpmath version: %s" % mpmath.__version__)
    print("Python version: %s" % sys.version)


def pytest_configure(config):
    config.addinivalue_line('markers', 'slow: marks tests as slow')


@pytest.fixture(autouse=True)
def reset_mp_globals():
    mpmath.mp.prec = sys.float_info.mant_dig
    mpmath.mp.pretty = False
    mpmath.mp.rounding = 'n'
    mpmath.mp.pretty_dps = "str"
    mpmath.iv.prec = mpmath.mp.prec
    mpmath.iv.pretty = False
