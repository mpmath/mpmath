import math

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
def fft(ctx, values, inverse=False, pad_to_power_of_two=False):
    r"""Computes the Fast Fourier Transform (or Inverse FFT) with optional zero-padding.

    If ``inverse=True``, computes the inverse transform and normalizes the
    result by the sequence length.
    If ``pad_to_power_of_two=True``, the input sequence is automatically zero-padded to the next power of 2 if its length is not already a power of 2.
    """
    orig_n = len(values)
    if orig_n == 0:
        return []

    is_power_of_two = (orig_n & (orig_n - 1)) == 0

    if not is_power_of_two:
        if not pad_to_power_of_two:
            raise ValueError(f"Length of input ({orig_n}) must be a power of 2.")

        padded_n = 1 << math.ceil(math.log2(orig_n))
    else:
        padded_n = orig_n

    converted_values = [ctx.convert(v) for v in values]

    if padded_n > orig_n:
        converted_values.extend([ctx.zero] * (padded_n - orig_n))

    result = _fft_recursive(ctx, converted_values, inverse)

    if inverse:
        scale = ctx.convert(padded_n)
        result = [val / scale for val in result]

    return result

@defun
def ifft(ctx, values, pad_to_power_of_two=False):
    r"""Computes the inverse discrete Fourier transform of a sequence with optional zero-padding.
    If ``pad_to_power_of_two=True``, the input sequence is automatically zero-padded to the next power of 2 if its length is not already a power of 2.
    """
    return fft(ctx, values, inverse=True, pad_to_power_of_two=pad_to_power_of_two)
