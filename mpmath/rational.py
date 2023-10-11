import fractions
import warnings


warnings.warn('mpmath.rational module is deprecated', DeprecationWarning)


def __getattr__(name):
    if name in ('mpq',):
        warnings.warn('mpmath.rational.mpq type is deprecated', DeprecationWarning)
        return globals()['_mpq']
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


class _mpq(fractions.Fraction):
    def __new__(cls, numerator=0, denominator=None):
        if isinstance(numerator, tuple) and len(numerator) > 1:
            return cls(numerator[0], numerator[1])
        return super().__new__(cls, numerator, denominator)

    @property
    def _mpq_(self):
        return self.as_integer_ratio()
