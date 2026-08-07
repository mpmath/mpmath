import math

from .calculus import defun

def _fft_iterative(ctx, values, inverse=False):
    """
    This function implements the Cooley-Tukey FFT algorithm iteratively. It computes the Fast Fourier Transform (or Inverse FFT) of a sequence of complex numbers.
    The input sequence must have a length that is a power of 2.

    https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
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

# TODO: Add support for non-power-of-two lengths using Bluestein's algorithm or similar methods.

@defun
def fft(ctx, values, pad_to_power_of_two=False):
    r"""Computes the Fast Fourier Transform (or Inverse FFT) with optional zero-padding.
    If ``pad_to_power_of_two=True``, the time domain sequence is automatically zero-padded at the end to the next power of 2 if its length is not already a power of 2.
    Raises NotImplementedError if the input sequence length is not a power of 2 and ``pad_to_power_of_two`` is False.
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
    If ``pad_to_power_of_two=True``, the frequency domain sequence is automatically center zero-padded with Nyquist-bin splitting to the next power of 2 if its length is not already a power of 2.
    The Nyquist-bin splitting is done by taking the value at the Nyquist frequency (if it exists) and splitting it evenly between the positive and negative Nyquist bins in the padded sequence.
    Raises NotImplementedError if the input sequence length is not a power of 2 and ``pad_to_power_of_two`` is False.
    The output is scaled by the length of the original input sequence to ensure that the inverse transform correctly recovers the original time-domain sequence.
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

    scale = ctx.convert(orig_n)
    result = [val / scale for val in result]

    return result
