def squarew(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> squarew(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        -1
    The average run time of the function using 200 dps (from 50 test runs) is 5.806399999999823e-05 seconds.

    """
    n = ctx.sinpi(2*t/T)
    if n == 0:
        return 0.0
    elif n > 0:
        return A
    else:
        # n < 0
        return -A


def squarew_floor(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> squarew_floor(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        -1.0
    The average run time of the function using 200 dps (from 50 test runs) is 2.081799999999967e-05 seconds.
    """
    return A*(4*ctx.floor(t/T) - 2*ctx.floor(2*t/T) + 1)


def squarew_floor_ex(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> squarew_floor_ex(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        -1.0
    The average run time of the function using 200 dps (from 50 test runs) is 1.2968000000001533e-05 seconds.

    """
    return A*((-1)**ctx.floor(2*t/T))


def trianglew(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> trianglew(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        -0.22222222222222232090871330001391470432281494140625
    The average run time of the function (from 50 test runs) using 200 dps is 0.00014768000000000281 seconds.

    """
    return (2*A/ctx.pi)*ctx.asin(ctx.sinpi(2*t/T))


def trianglew_saw(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> trianglew_saw(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        -0.22222222222222232090871330001391470432281494140625
    The average run time of the function (from 50 test runs) using 200 dps is 0.0001810279999999992 seconds.

    """
    n = t*ctx.pi/T
    if n == ctx.pi/4 or (n - ctx.pi/4) % ctx.pi == 0:
        return A*(abs((4/ctx.pi)*(ctx.pi/2)) - 1)
    else:
        return A*(abs((4/ctx.pi)*ctx.atan(ctx.cot((t*ctx.pi/T) - ctx.pi/4))) - 1)


def sawtoothw(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> sawtoothw(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        0.111111111111111160454356650006957352161407470703125
    The average run time of the function (from 50 test runs) using 200 dps is 0.0001736160000000009 seconds.

    """
    n = t*ctx.pi/T
    if n % ctx.pi == 0:
        return -(2*A/ctx.pi)*(ctx.pi/2)
    elif (n - ctx.pi/2) % ctx.pi == 0:
        return 0.0
    else:
        return -(2*A/ctx.pi)*ctx.atan(ctx.cot(t*ctx.pi/T))


def sawtoothw_mod(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> sawtoothw_mod(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        0.111111111111111160454356650006957352161407470703125
    The average run time of the function (from 50 test runs) using 200 dps is 1.6823999999998896e-05 seconds.

    """
    return A*((2/T)*ctx.fmod(t, T) - 1)


def sawtoothw_floor(ctx, t, A, T):
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

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> sawtoothw_floor(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1,2)
        0.111111111111111160454356650006957352161407470703125
    The average run time of the function (from 50 test runs) using 200 dps is 1.763999999999988e-05 seconds.

    """
    return A*(2*((t/T) - ctx.floor(t/T)) - 1)


def sigmoidw(ctx, t, A):
    """
    Computes the sigmoid wave function using the definition:

    .. math::
        x(t) = \frac{A}{1 + e^{-x}} - A

    where A is the amplitude.

    Testing run time with 200 dps ::
        >>> from mpmath import *
        >>> mp.dps = 200; mp.pretty = True
        >>> sigmoidw(1.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111,1)
        0.50467239772185676840241452244896871205195656825700120504353593216961978482423997202634746804585666419751329795595729358385128816779304161122766839674446523241806725207242033774310388431511520878680774
    The average run time of the function (from 50 test runs) using 200 dps is 0.00019384400000000023 seconds.

    """
    return 2*A/(1 + ctx.exp(-t)) - A
