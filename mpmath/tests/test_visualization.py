"""
Limited tests of the visualization module. Right now it just makes
sure that passing custom Axes works.

"""

import pytest

from mpmath import fp, mp


def test_axes():
    try:
        import matplotlib
        version = matplotlib.__version__.split("-")[0]
        version = version.split(".")[:2]
        if [int(_) for _ in version] < [0,99]:
            raise ImportError
        import pylab
    except ImportError:
        pytest.skip("\nSkipping test (pylab not available or too old version)\n")
    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for ctx in [mp, fp]:
        ctx.plot(lambda x: x**2, [0, 3], axes=axes)
        assert axes.get_xlabel() == 'x'
        assert axes.get_ylabel() == 'f(x)'

    fig = pylab.figure()
    axes = fig.add_subplot(111)
    for ctx in [mp, fp]:
        ctx.cplot(lambda z: z, [-2, 2], [-10, 10], axes=axes)
    assert axes.get_xlabel() == 'Re(z)'
    assert axes.get_ylabel() == 'Im(z)'


def test_issue_1007():
    # plot(), cplot() and splot() must not leave a stale figure open
    # when the user-supplied function raises an unexpected exception;
    # otherwise that blank figure lingers and is shown on the next call.
    try:
        import matplotlib
        version = matplotlib.__version__.split("-")[0]
        version = version.split(".")[:2]
        if [int(_) for _ in version] < [0,99]:
            raise ImportError
        import pylab
    except ImportError:
        pytest.skip("\nSkipping test (pylab not available or too old version)\n")

    class Boom(Exception):
        pass

    def bad(*args):
        # An error that is not in plot_ignore, so it propagates out of
        # plot()/cplot()/splot() instead of being silently skipped.
        raise Boom

    for ctx in [mp, fp]:
        pylab.close("all")
        pytest.raises(Boom, lambda: ctx.plot(bad, [0, 2]))
        assert pylab.get_fignums() == []

        pylab.close("all")
        pytest.raises(Boom, lambda: ctx.cplot(bad, [-2, 2], [-2, 2]))
        assert pylab.get_fignums() == []

        pylab.close("all")
        pytest.raises(Boom, lambda: ctx.splot(bad, [-1, 1], [-1, 1]))
        assert pylab.get_fignums() == []
