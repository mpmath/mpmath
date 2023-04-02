"""
Modules checks, in particular some hacks for Issue #657
"""

from pytest import raises

import mpmath


def test_erroneous_module_setting():
    with raises(AttributeError):
        mpmath.dps = 64
    with raises(AttributeError):
        mpmath.prec = 192
    with raises(AttributeError):
        mpmath.pretty = True
    with raises(AttributeError):
        mpmath.trap_complex = True
