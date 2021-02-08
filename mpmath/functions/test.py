import mpmath

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


def trianglew(t: float, A: float, T: float) -> float:
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

mpmath.mp.dps = 200; mpmath.mp.pretty = True

for i in range(-1000,-995):

    temp = mpmath.fdiv(i,100)
    print(temp)
    answer = trianglew(temp,1,1)
    print(answer)
    answer = trianglew_saw(temp, 1, 1)
    print(answer)
    print(" ")