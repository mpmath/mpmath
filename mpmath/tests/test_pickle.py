import pickle

import pytest

from mpmath import matrix, mpc, mpf, mpi, sin


@pytest.mark.parametrize('protocol', range(pickle.HIGHEST_PROTOCOL + 1))
@pytest.mark.parametrize('obj', [mpf('0.5'), mpc('0.5','0.2'), mpi(10, 30),
                                 matrix([1, sin(1)]), matrix([[1, 2], [3, 4]])])
def test_pickle(obj, protocol):
    assert obj == pickle.loads(pickle.dumps(obj, protocol))
