import pytest

from mpmath import fp, zetazero


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
