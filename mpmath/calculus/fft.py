from .calculus import defun

def _fft_recursive(ctx, values, inverse=False):
    n = len(values)

    if n == 1:
        return values

    sign = 1 if inverse else -1

    even = _fft_recursive(ctx, values[0::2], inverse)
    odd = _fft_recursive(ctx, values[1::2], inverse)
    root = ctx.exp(sign * 2 * ctx.pi * ctx.j / n)
    factor = ctx.one
    half = n // 2
    transformed = [ctx.zero] * n
    for k in range(half):
        term = factor * odd[k]
        transformed[k] = even[k] + term
        transformed[k + half] = even[k] - term
        factor *= root
    return transformed

@defun
def fft(ctx, values, inverse=False):
    r"""Computes the discrete Fourier transform of a sequence.

    If ``inverse=True``, computes the inverse transform and normalizes the
    result by the sequence length.
    """
    n = len(values)
    if n == 0:
        return []
    if n & (n - 1) != 0:
        raise ValueError("Length of input must be a power of 2.")

    converted_values = [ctx.convert(v) for v in values]
    transformed = _fft_recursive(ctx, converted_values, inverse=inverse)
    if inverse and transformed:
        n = len(transformed)
        transformed = [value / n for value in transformed]
    return transformed

@defun
def ifft(ctx, values):
    r"""Computes the inverse discrete Fourier transform of a sequence."""
    return fft(ctx, values, inverse=True)
