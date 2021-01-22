# order of square wave function in terms of increasing run time:
# square_wave_floor_ex, square_wave, square_wave_floor


def squarew(ctx,t, A, T):
    """
    Computes the square wave function using the definition:

    ..math::
        x(t) = Asgn(\sin(\frac{2 \pi t}{T}))

    where sgn is the sign function and x is 0 at the discontinuities.

    **Examples**

    Square wave with period = 2, amplitude = 1 ::

        >>> squarew(0,1,2)
        0.0
        >>> squarew(0.5,1,2)
        1
        >>> squarew(1,1,2)
        0.0
        >>> squarew(1.5,1,2)
        -1
        >>> squarew(2,1,2)
        0.0

    """
    n = ctx.sinpi(2*t/T)
    if n == 0:
        return 0.0
    elif n > 0:
        return A
    else:
        # n < 0
        return -A


def squarew_floor(ctx,t, A, T):
    """
    Computes the square wave function using the definition:

    .. math::
        x(t) = A(4\left\lfloor{\frac{t}{T}}\right\rfloor - 2\left\lfloor{\frac{2t}{T}}\right\rfloor + 1)

    where T is the period of the wave and A is the amplitude.

    **Examples**

    Square wave with period = 2, amplitude = 1 ::

        >>> squarew_floor(0,1,2)
        1.0
        >>> squarew_floor(0.5,1,2)
        1.0
        >>> squarew_floor(1,1,2)
        -1.0
        >>> squarew_floor(1.5,1,2)
        -1.0
        >>> squarew_floor(2,1,2)
        1.0
    """
    return A*(4*ctx.floor(t/T) - 2*ctx.floor(2*t/T) + 1)


def squarew_floorex(ctx,t, A, T):
    """
    Computes the square wave function using the definition:

    .. math::
        x(t) = A(-1)^{\left\lfloor{\frac{2t}{T}}\right\rfloor}

    where T is the period of the wave and A is the amplitude.

    **Examples**

    Square wave with period = 2, amplitude = 1 ::

        >>> squarew_floor_ex(0,1,2)
        1.0
        >>> squarew_floor_ex(0.5,1,2)
        1.0
        >>> squarew_floor_ex(1,1,2)
        -1.0
        >>> squarew_floor_ex(1.5,1,2)
        -1.0
        >>> squarew_floor_ex(2,1,2)
        1.0

    """
    return A*((-1)**ctx.floor(2*t/T))

# order of triangle wave function in terms of increasing run time:
# triangle_wave, triangle_wave_saw

def trianglew(ctx,t, A, T):
    """
    Computes the triangle wave function using the definition:

    .. math::
        x(t) = \frac{2A}{\pi}\arcsin(\sin(\frac{2\pi t}{T}))

    where T is the period of the wave and A is the amplitude.

    **Examples**

    Triangle wave with period = 2, amplitude = 1 ::

        >>> trianglew(0,1,2)
        0.0
        >>> trianglew(0.25,1,2)
        0.5
        >>> trianglew(0.5,1,2)
        1.0
        >>> trianglew(1,1,2)
        0.0
        >>> trianglew(1.5,1,2)
        -1.0
        >>> trianglew(2,1,2)
        0.0
    """
    return (2*A/ctx.pi)*ctx.asin(ctx.sinpi(2*t/T))


def trianglew_saw(ctx,t, A, T):
    """
    Computes the triangle wave function as the absolute function of a sawtooth wave:

    .. math::
        x(t) = A(|(\frac{4}{\pi})\arctan(\cot(\frac{\pi t}{T} - \frac{\pi}{4}))| - 1)

    where T is the period of the wave, and A is the amplitude.

    **Examples**

    Triangle wave with period = 2, amplitude = 1 ::

        >>> trianglew_saw(0,1,2)
        1.0
        >>> trianglew_saw(0.25,1,2)
        0.5
        >>> trianglew_saw(0.5,1,2)
        1.0
        >>> trianglew_saw(1,1,2)
        0.0
        >>> trianglew_saw(1.5,1,2)
        -1.0
        >>> trianglew_saw(2,1,2)
        1.0
    """
    n = t*ctx.pi/T
    if n == ctx.pi/4 or (n - ctx.pi/4) % ctx.pi == 0:
        return A*(abs((4/ctx.pi)*(ctx.pi/2)) - 1)
    else:
        return A*(abs((4/ctx.pi)*ctx.atan(ctx.cot((t*ctx.pi/T) - ctx.pi/4))) - 1)


def sawtoothw(t, A, T):
    """
    Computes the sawtooth wave function using the definition:

    .. math::
        x(t) = -(\frac{2A}{\pi})\arctan(\cot(\frac{\pi t}{T}))

    where T is the period and A is the amplitude.

    **Examples**

    Sawtooth wave with period = 2, amplitude = 1 ::

        >>> sawtoothw(0,1,2)
        -1.0
        >>> sawtoothw(0.25,1,2)
        -0.75
        >>> sawtoothw(0.5,1,2)
        -0.5
        >>> sawtoothw(1,1,2)
        0.0
        >>> sawtoothw(1.5,1,2)
        0.5
        >>> sawtoothw(2,1,2)
        -1.0

    """
    n = t*ctx.pi/T
    if n % ctx.pi == 0:
        return -(2*A/ctx.pi)*(ctx.pi/2)
    elif (n - ctx.pi/2) % ctx.pi == 0:
        return 0.0
    else:
        return -(2*A/ctx.pi)*ctx.atan(ctx.cot(t*ctx.pi/T))

# order of sawtooth wave functions in increasing run time:
# sawtooth_wave_mod, sawtooth_wave_floor, sawtooth_wave

def sawtoothw_mod(t, A, T):
    """
    Computes the sawtooth wave function using the definition:

    .. math::
        x(t) = A(\frac{2}{T}(t \bmod T) - 1)

    where T is the period and A is the amplitude.

    **Examples**

    Sawtooth wave with period = 2, amplitude = 1 ::

        >>> sawtoothw_mod(0,1,2)
        -1.0
        >>> sawtoothw_mod(0.25,1,2)
        -0.75
        >>> sawtoothw_mod(0.5,1,2)
        -0.5
        >>> sawtoothw_mod(1,1,2)
        0.0
        >>> sawtoothw_mod(1.5,1,2)
        0.5
        >>> sawtoothw_mod(2,1,2)
        -1.0

    """
    return A*((2/T)*ctx.fmod(t, T) - 1)


def sawtoothw_floor(t, A, T):
    """
    Computes the sawtooth wave function using the definition:

    .. math::
        x(t) = A(2(\frac{t}{T} - \left\lfloor{\frac{t}{T}}\right\rfloor) - 1)

    where T is the period and A is the amplitude.

    **Examples**

    Sawtooth wave with period = 2, amplitude = 1 ::

        >>> sawtoothw_floor(0,1,2)
        -1.0
        >>> sawtoothw_floor(0.25,1,2)
        -0.75
        >>> sawtoothw_floor(0.5,1,2)
        -0.5
        >>> sawtoothw_floor(1,1,2)
        0.0
        >>> sawtoothw_floor(1.5,1,2)
        0.5
        >>> sawtoothw_floor(2,1,2)
        -1.0

    """
    return A*(2*((t/T) - ctx.floor(t/T)) - 1)

def sigmoidw(t, A):
    """
    Computes the sigmoid wave function using the definition:

    .. math::
        x(t) = \frac{A}{1 + e^{-x}} - A

    where A is the amplitude.
    """
    return 2A/(1 + exp(-t)) - A
