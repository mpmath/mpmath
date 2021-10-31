
def squarew(ctx, t, amplitude=1, period=1):
    r"""
    Computes the square wave function using the definition:

    .. math::
        x(t) = A(-1)^{\left\lfloor{\frac{2t}{P}}\right\rfloor}

    where `P` is the period of the wave and `A` is the amplitude.

    **Examples**

    Square wave with period = 2, amplitude = 1 ::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> squarew(0,1,2)
        1.0
        >>> squarew(0.5,1,2)
        1.0
        >>> squarew(1,1,2)
        -1.0
        >>> squarew(1.5,1,2)
        -1.0
        >>> squarew(2,1,2)
        1.0
    """
    P = period
    A = amplitude
    return A*((-1)**ctx.floor(2*t/P))


def trianglew(ctx, t, amplitude=1, period=1):
    r"""
    Computes the triangle wave function using the definition:

    .. math::
        x(t) = 2A\left(\frac{1}{2}-\left|1-2frac\left(\frac{x}{P}+\frac{1}{4}\right)\right|\right)

    where :math:`frac\left(\frac{t}{T}\right) = \frac{t}{T}-\left\lfloor{\frac{t}{T}}\right\rfloor`
    , `P` is the period of the wave, and `A` is the amplitude.

    **Examples**

    Triangle wave with period = 2, amplitude = 1 ::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
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
    A = amplitude
    P = period

    return 2*A*(0.5 - ctx.fabs(1 - 2*ctx.frac(t/P + 0.25)))


def sawtoothw(ctx, t, amplitude=1, period=1):
    r"""
    Computes the sawtooth wave function using the definition:

    .. math::
        x(t) = Afrac\left(\frac{t}{T}\right)

    where :math:`frac\left(\frac{t}{T}\right) = \frac{t}{T}-\left\lfloor{\frac{t}{T}}\right\rfloor`
    , `P` is the period of the wave, and `A` is the amplitude.

    **Examples**

    Sawtooth wave with period = 2, amplitude = 1 ::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> sawtoothw(0,1,2)
        0.0
        >>> sawtoothw(0.5,1,2)
        0.25
        >>> sawtoothw(1,1,2)
        0.5
        >>> sawtoothw(1.5,1,2)
        0.75
        >>> sawtoothw(2,1,2)
        0.0
    """
    A = amplitude
    P = period
    return A*ctx.frac(t/P)


def unit_triangle(t, amplitude=1):
    r"""
    Computes the unit triangle using the definition:

    .. math::
        x(t) = A(-\left| t \right| + 1)

    where `A` is the amplitude.

    **Examples**

    Unit triangle with amplitude = 1 ::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> unit_triangle(-1,1)
        0
        >>> unit_triangle(-0.5,1)
        0.5
        >>> unit_triangle(0,1)
        1
        >>> unit_triangle(0.5,1)
        0.5
        >>> unit_triangle(1,1)
        0
    """
    A = amplitude
    if t <= -1 or t >= 1:
        return ctx.zero
    return A*(-abs(t) + 1)


def sigmoid(ctx, t, amplitude=1):
    r"""
    Computes the sigmoid function using the definition:

    .. math::
        x(t) = \frac{A}{1 + e^{-t}}

    where `A` is the amplitude.

    **Examples**

    Sigmoid function with amplitude = 1 ::

        >>> from mpmath import *
        >>> mp.dps = 25; mp.pretty = True
        >>> sigmoid(-1,1)
        0.2689414213699951207488408
        >>> sigmoid(-0.5,1)
        0.3775406687981454353610994
        >>> sigmoid(0,1)
        0.5
        >>> sigmoid(0.5,1)
        0.6224593312018545646389006
        >>> sigmoid(1,1)
        0.7310585786300048792511592

    """
    A = amplitude
    return A / (1 + ctx.exp(-t))
