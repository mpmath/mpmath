import math

from .calculus import defun

def _fft_iterative(ctx, values, inverse=False):
    """
    This function implements the Cooley-Tukey FFT algorithm iteratively. It computes the Fast Fourier Transform (or Inverse FFT) of a sequence of complex numbers.
    The input sequence must have a length that is a power of 2. If the length is not a power of 2, the function raises a NotImplementedError.

    https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Data_reordering,_bit_reversal,_and_in-place_algorithms
    """
    n = len(values)
    if n <= 1:
        return values

    # Bit-Reversal Permutation
    transformed = [ctx.zero] * n
    num_bits = n.bit_length() - 1
    for i in range(n):
        rev = 0
        for b in range(num_bits):
            if (i >> b) & 1:
                rev |= 1 << (num_bits - 1 - b)
        transformed[rev] = values[i]

    sign = 1 if inverse else -1
    two_pi_j = 2 * ctx.pi * ctx.j

    length = 2
    while length <= n:
        half = length // 2
        w_len = ctx.exp(sign * two_pi_j / length)

        for i in range(0, n, length):
            w = ctx.one
            for j in range(half):
                u = transformed[i + j]
                v = transformed[i + j + half] * w

                transformed[i + j] = u + v
                transformed[i + j + half] = u - v

                w *= w_len

        length <<= 1

    return transformed

@defun
def fft(ctx, values, pad_to_power_of_two=False):
    r"""Computes the Fast Fourier Transform (or Inverse FFT) with optional zero-padding.
    If ``pad_to_power_of_two=True``, the time domain sequence is automatically zero-padded at the end to the next power of 2 if its length is not already a power of 2.
    """
    orig_n = len(values)
    if orig_n == 0:
        return []

    is_power_of_two = (orig_n & (orig_n - 1)) == 0

    if not is_power_of_two:
        if not pad_to_power_of_two:
            raise NotImplementedError("FFT is only implemented for lengths that are powers of 2. Use pad_to_power_of_two=True to automatically pad the input sequence.")

        padded_n = 1 << math.ceil(math.log2(orig_n))
    else:
        padded_n = orig_n

    converted_values = [ctx.convert(v) for v in values]
    pad_length = padded_n - orig_n

    if padded_n > orig_n:
        converted_values.extend([ctx.zero] * pad_length)

    result = _fft_iterative(ctx, converted_values, False)

    return result

@defun
def invfft(ctx, values, pad_to_power_of_two=False):
    r"""Computes the inverse discrete Fourier transform of a sequence with optional zero-padding.
    If ``pad_to_power_of_two=True``, the frequency domain sequence is automatically zero-padded with Nyquist-bin splitting to the next power of 2 if its length is not already a power of 2.
    """

    orig_n = len(values)
    if orig_n == 0:
        return []

    is_power_of_two = (orig_n & (orig_n - 1)) == 0

    if not is_power_of_two:
        if not pad_to_power_of_two:
            raise NotImplementedError("Inverse FFT is only implemented for lengths that are powers of 2. Use pad_to_power_of_two=True to automatically pad the input sequence.")

        padded_n = 1 << math.ceil(math.log2(orig_n))
    else:
        padded_n = orig_n

    converted_values = [ctx.convert(v) for v in values]
    pad_length = padded_n - orig_n

    half = orig_n >> 1

    if pad_length > 0:
        if orig_n % 2 == 0:
            zeros = [ctx.zero] * (pad_length - 1)
            nyquist_half = converted_values[half] / 2
            converted_values = (
                converted_values[:half]
                + [nyquist_half]
                + zeros
                + [nyquist_half]
                + converted_values[half + 1:]
                )
        else:
            zeros = [ctx.zero] * pad_length
            converted_values = (
                converted_values[:half + 1]
                + zeros
                + converted_values[half + 1:]
                )

    result = _fft_iterative(ctx, converted_values, True)

    scale = ctx.convert(padded_n)
    result = [val / scale for val in result]

    return result
