import pytest

from mpmath import mp
from mpmath.libmp import from_int, round_ceiling, round_floor, round_nearest
from mpmath.libmp.gammazeta import mpf_psi, mpf_psi0


@pytest.mark.parametrize('prec', [15, 53, 113])
@pytest.mark.parametrize('value', [1, 3, '0.75', complex(1.25, 0.5)])
def test_psi_zero_is_rounded_to_context_precision(prec, value):
    ctx = mp.clone()
    ctx.prec = prec
    z = ctx.mpc(value) if isinstance(value, complex) else ctx.mpf(value)
    result = ctx.psi(0, z)
    assert result == +result
    if isinstance(result, ctx.mpc):
        assert result.real._mpf_[3] <= prec
        assert result.imag._mpf_[3] <= prec
    else:
        assert result._mpf_[3] <= prec


@pytest.mark.parametrize('value', [1, 3, '0.75'])
def test_psi_zero_matches_digamma_rounding(value):
    ctx = mp.clone()
    ctx.prec = 53
    z = ctx.mpf(value)
    assert ctx.psi(0, z) == ctx.digamma(z)


@pytest.mark.parametrize('rounding', [round_ceiling, round_floor, round_nearest])
def test_real_psi_zero_preserves_requested_rounding(rounding):
    value = from_int(3)
    assert mpf_psi(0, value, 53, rounding) == mpf_psi0(value, 53, rounding)