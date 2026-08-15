from .calculus import defun

def _fft_cooley_tuckey(ctx, values, inverse=False):
    """
    This function implements the Radix-2 Cooley-Tukey FFT algorithm iteratively. It computes the Fast Fourier Transform (or Inverse FFT) of a sequence of complex numbers.
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
def fft(ctx, values):
    r"""Computes the Discrete Fourier Transform (DFT) of a sequence using the Radix-2 Cooley-Tukey FFT algorithm.
    Raises NotImplementedError if the input sequence length is not a power of 2.
    """
    n = len(values)
    if n == 0:
        return []

    is_power_of_two = (n & (n - 1)) == 0

    if not is_power_of_two:
        raise NotImplementedError(f"FFT is only implemented for lengths that are powers of 2, got length: {n}")

    converted_values = [ctx.convert(v) for v in values]

    result = _fft_cooley_tuckey(ctx, converted_values, False)

    return result

@defun
def invfft(ctx, values):
    r"""Computes the inverse Discrete Fourier Transform (IDFT) of a sequence using the Radix-2 Cooley-Tukey FFT algorithm.
    Raises NotImplementedError if the input sequence length is not a power of 2.
    """

    n = len(values)
    if n == 0:
        return []

    is_power_of_two = (n & (n - 1)) == 0

    if not is_power_of_two:
        raise NotImplementedError(f"Inverse FFT is only implemented for lengths that are powers of 2, got length: {n}")

    converted_values = [ctx.convert(v) for v in values]

    result = _fft_cooley_tuckey(ctx, converted_values, True)

    result = [val / n for val in result]

    return result
