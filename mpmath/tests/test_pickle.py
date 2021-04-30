import os
import tempfile
import pickle

from mpmath import *

def pickler(obj, protocol=0):
    fn = tempfile.mktemp()

    with open(fn, 'wb') as f:
        pickle.dump(obj, f, protocol=protocol)

    with open(fn, 'rb') as f:
        obj2 = pickle.load(f)

    os.remove(fn)

    return obj2

def test_pickle():

    obj = mpf('0.5')

    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        assert obj == pickler(obj, protocol)

    obj = mpc('0.5','0.2')
    for protocol in range(pickle.HIGHEST_PROTOCOL + 1):
        assert obj == pickler(obj, protocol)
