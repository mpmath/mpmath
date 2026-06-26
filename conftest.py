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

    if "no:hypothesispytest" not in config.getoption("-p"):
        from hypothesis import settings

        default = settings.get_profile("default")
        settings.register_profile("default",
                                  settings(default, max_examples=1000))
        ci = settings.get_profile("ci")
        settings.register_profile("ci", settings(ci, max_examples=10000))


@pytest.fixture(autouse=True)
def reset_mp_globals():
    mpmath.mp.prec = sys.float_info.mant_dig
    mpmath.mp.pretty = False
    mpmath.mp.rounding = 'n'
    mpmath.mp.pretty_dps = "str"
    mpmath.mp._legacy = True
    mpmath.iv.prec = mpmath.mp.prec
    mpmath.iv.pretty = False
