from mpmath import mp
import mpmath
import numpy as np

mpmath.mp.dps = 200

def squarew(t: float, A: float, T: float) -> float:
    n = mpmath.sinpi(2*t/T)
    if n == 0:
        return 0.0
    elif n > 0:
        return A
    else:
        # n < 0
        return -A


def squarew_floor(t: float, A: float, T: float) -> float:
    return A*(4*mpmath.floor(t/T) - 2*mpmath.floor(2*t/T) + 1)


def squarew_floor_ex(t: float, A: float, T: float) -> float:
    return A*((-1)**mpmath.floor(2*t/T))


def trianglew(
    t: float,
    A: float = 1,
    T: float = 1,
) -> float:
    return (2*A/mpmath.pi)*mpmath.asin(mpmath.sinpi(2*t/T))

def trianglew_saw(t: float, A: float, T: float) -> float:
    n = t*mpmath.pi/T
    if n == mpmath.pi/4 or (n - mpmath.pi/4) % mpmath.pi == 0:
        return A*(abs((4/mpmath.pi)*(mpmath.pi/2)) - 1)
    else:
        return A*(abs((4/mpmath.pi)*mpmath.atan(mpmath.cot((t*mpmath.pi/T) - mpmath.pi/4))) - 1)


def sawtoothw(t: float, A: float, T: float) -> float:
    n = t*mpmath.pi/T
    if n % mpmath.pi == 0:
        return -(2*A/mpmath.pi)*(mpmath.pi/2)
    elif (n - mpmath.pi/2) % mpmath.pi == 0:
        return 0.0
    else:
        return -(2*A/mpmath.pi)*mpmath.atan(mpmath.cot(t*mpmath.pi/T))


def sawtoothw_mod(t: float, A: float, T: float) -> float:
    return A*((2/T)*mpmath.fmod(t, T) - 1)


def sawtoothw_floor(t: float, A: float, T: float) -> float:
    return A*(2*((t/T) - mpmath.floor(t/T)) - 1)


def unit_triangle(t: float, A: float) -> float:
    if t <= -1 or t >= 1:
        return 0
    else:
        return A*(-abs(t) + 1)


def sigmoidw(t: float, A: float) -> float:
    return A/(1 + mpmath.exp(-t))

def triangle_wave(x, frange=(-1, 1), period=1):
    A = frange[1] - frange[0]
    bottom = frange[0]
    P = period

    return A*(1-2*mp.fabs(0.5-mp.frac(x/P+0.25))) + bottom

N = -10
step = mp.fdiv(1, 100)

while N <= -9.9:
    x = N
    print(float(x))
    answer = sigmoidw(x, 1)
    
    calib_row = calibration_data[calibration_data[:, 0].astype(np.float) == float(x)].copy()
    calib_value = mp.mpf(calib_row[0,5])
    
    print(mp.fsub(answer, calib_value))
    print("")
    
    N += step