import pytest


def pytest_configure(config):
    config.addinivalue_line('markers', 'slow: marks tests as slow')


@pytest.fixture(autouse=True)
def reset_mp_globals():
    from mpmath import mp, iv
    mp.dps = 15
    mp.pretty = False
    iv.dps = 15
    iv.pretty = False
