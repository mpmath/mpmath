import pytest

from mpmath import fp, mp, mpf, zeta, zetazero


@pytest.mark.parametrize("n,v",
    [(399999999, 156762524.6750591511),
     (241389216, 97490234.2276711795),
     (526196239, 202950727.691229534),
     (542964976, 209039046.578535272),
     (1048449112, 388858885.231056486),
     (1048449113, 388858885.384337406),
     (1048449114, 388858886.002285122),
     (1048449115, 388858886.00239369),
     (1048449116, 388858886.690745053),
     (3570918901, 1239587702.54745031),
     (3570918902, 1239587702.54752387),
     # issue 1147, see
     # https://www.lmfdb.org/zeros/zeta/?limit=100&N=325890640
     # and https://www.lmfdb.org/zeros/zeta/?limit=100&N=325890640
     (325890640, 129273228.66142665),
     (325890641, 129273228.76005181),
     (325890642, 129273228.79754069),
     (357738764, 141125096.01260684),
     (357738765, 141125096.18511831),
     (357738766, 141125096.28064566),
     # Huge zeros (this may take hours):
#    (8637740722917, 2124447368584.39296466152),
#    (8637740722918, 2124447368584.39298170604),
     ])
def test_zetazero(n, v):
    assert zetazero(n).ae(complex(0.5,v))

def test_zeta_param(capsys):
    fp.zeta(0.5+100j, method="riemann-siegel", verbose=True)
    captured = capsys.readouterr()
    assert "Attempting to use the Riemann-Siegel algorithm" in captured.out


def test_hurwitz_zeta_tiny_values():
    # Issue #1198: the Euler-Maclaurin tail in _hurwitz_em stopped on an
    # absolute 2**-prec tolerance, so values far below 2**-prec were only
    # accurate to a few digits (2.8e-9 relative error for this case at 30 dps).
    # Inputs are exact binary floats; the reference is direct summation of
    # (q+k)**-s plus the integral tail, at 60 dps.
    s = mpf(43.01303985093054)
    q = mpf(305.4243263250737)
    mp.dps = 30
    ref = mpf('1.018131113260128483346694e-106')
    assert abs(zeta(s, q) - ref) / ref < mpf(10)**-24
    mp.dps = 15
    assert abs(zeta(s, q) - ref) / ref < mpf(10)**-13
