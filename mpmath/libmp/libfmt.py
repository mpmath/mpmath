from .libmpf import (to_digits_exp, stddigits, round_nearest, fzero, fnzero,
                     finf, fninf, fnan)
blog2 = 3.3219280948873626


def round_digits(digits, dps, base, rounding=round_nearest):
    '''
    Returns the rounded digits, and the number of places the decimal point was
    shifted.
    '''

    exponent = 0

    if rounding == round_nearest:
        rnd_digs = stddigits[(base//2 + base % 2):base]
    else:
        raise ValueError('Rounding: {} not implemented yet.')

    round_up = False
    if len(digits) > dps:
        # The first digit after dps is a 5.
        if digits[dps] == rnd_digs[0]:
            for i in range(dps+1, len(digits)):
                if digits[i] != '0':
                    round_up = True
                    break
            if digits[dps-1] in stddigits[1:base:2]:
                round_up = True

    if not round_up:
        digits = digits[:dps] + stddigits[int(digits[dps], base) - 1]

    # Rounding up kills some instances of "...99999"
    if len(digits) > dps and digits[dps] in rnd_digs:
        digits = digits[:dps]
        i = dps - 1
        dig = stddigits[base-1]
        while i >= 0 and digits[i] == dig:
            i -= 1
        if i >= 0:
            digits = digits[:i] + stddigits[int(digits[i], base) + 1] + \
                '0' * (dps - i - 1)
        else:
            digits = '1' + '0' * (dps - 1)
            exponent += 1
    else:
        digits = digits[:dps]

    return digits, exponent


def calc_padding(nchars, width, align):
    '''
    Computes the left and right padding required to fill the required width,
    according to how the string will be aligned.
    '''
    ntotal = max(nchars, width)

    if align == '<':
        lpad = 0
        rpad = ntotal - nchars
    elif align in ('>', '='):
        lpad = ntotal - nchars
        rpad = 0
    elif align == '^':
        lpad = (ntotal - nchars)//2
        rpad = ntotal - nchars - lpad
    else:
        lpad = 0
        rpad = 0

    return (lpad, rpad)


def read_format_spec(format_spec):
    '''
    Reads the format spec into a dictionary.
    This is more or less copied from the CPython implementation for regular
    floats.
    '''

    pos = 0

    format_dict = {
        'fill_char': ' ',
        'align': '>',
        'sign': '-',
        'no_neg_0': False,
        'alternate': False,
        'thousands_separators': '',
        'width': -1,
        'precision': 6,
        'rounding': 'N',
        'type': 'f'
        }

    # Check for alignment and fill characters
    if pos+1 < len(format_spec):
        if format_spec[pos+1] in ('<', '^', '>', '='):
            format_dict['fill_char'] = format_spec[pos]
            format_dict['align'] = format_spec[pos+1]
            pos += 2

    # Check for only alignment characters
    if pos < len(format_spec):
        if format_spec[pos] in ('<', '^', '>', '='):
            format_dict['fill_char'] = ' '
            format_dict['align'] = format_spec[pos]
            pos += 1

    # Check for sign character
    if pos < len(format_spec):
        if format_spec[pos] in ('+', '-', ' '):
            format_dict['sign'] = format_spec[pos]
            pos += 1

    # Check if negating zero
    if pos < len(format_spec):
        if format_spec[pos] == 'z':
            format_dict['no_neg_0'] = True
            pos += 1

    # Check if alternate
    if pos < len(format_spec):
        if format_spec[pos] == '#':
            format_dict['alternate'] = True
            pos += 1

    # Check for width
    if pos < len(format_spec):
        dig_str = ''
        while format_spec[pos].isdigit():
            dig_str += format_spec[pos]
            pos += 1

        if len(dig_str) > 0:
            format_dict['width'] = int(dig_str)

    # Check for thousands separator
    if pos < len(format_spec):
        if format_spec[pos] in (',', '_'):
            format_dict['thousands_separators'] = format_spec[pos]
            pos += 1

    # Check for precision
    if pos < len(format_spec):
        if format_spec[pos] == '.':
            pos += 1
            dig_str = ''
            while format_spec[pos].isdigit():
                dig_str += format_spec[pos]
                pos += 1
            if len(dig_str) == 0:
                raise ValueError('Precision not specified')
            else:
                format_dict['precision'] = int(dig_str)

    # Check for rounding type
    if pos < len(format_spec):
        if format_spec[pos] in ('U', 'D', 'Y', 'Z', 'N'):
            format_dict['rounding'] = format_spec[pos]
            pos += 1

    # Check for format type
    if pos < len(format_spec):
        if format_spec[pos] in ('f', 'F', 'g', 'G', 'e', 'E'):
            format_dict['type'] = format_spec[pos]
        else:
            raise ValueError(
                    "Format type '{}' not recognized.".format(format_spec[pos])
                    )
        pos += 1

    # -------------------------------------------------------------------------
    # Now process information

    # Compute actual width
    return format_dict


def format_fixed(s,
                 precision=6,
                 strip_zeros=False,
                 strip_last_zero=False,
                 thousands_separators='',
                 sign_spec='-',
                 sep_range=3,
                 base=10,
                 alternate=False):
    '''
    Format a number into fixed point.
    Returns the sign character, and the string that represents the number in
    the correct format.
    Does not perform padding or aligning
    '''

    # First, get the exponent to know how many digits we will need
    _, _, exponent = to_digits_exp(s, 1, base)

    # Now that we have an estimate, compute the correct digits
    # (we do this because the previous computation could yield the wrong
    # exponent by +- 1)
    sign, digits, exponent = to_digits_exp(
            s, max(precision+exponent+4, int(s[3]/blog2)), base)
    dps = precision + exponent + 1

    # Hack: if the digits are all 9s, then we will lose one dps when rounding
    # up.
    if all(dig == stddigits[base-1] for dig in digits[:dps+1]):
        dps += 1

    if sign != '-' and sign_spec != '-':
        sign = sign_spec

    # The number we want to print is lower in magnitude that the requested
    # precision. We should only print 0s.
    if exponent + precision < -1:
        if precision > 0:
            digits = '0.' + precision*'0'
        else:
            digits = '0'
    else:
        # Rounding up kills some instances of "...99999"
        digits, exp_add = round_digits(digits, dps, base)
        exponent += exp_add

        # Here we prepend the corresponding 0s to the digits string, according
        # to the value of exponent
        if exponent < 0:
            digits = ("0"*(-exponent)) + digits
            split = 1
        else:
            split = exponent + 1
            if split > dps:
                digits += "0"*(split-dps)
        exponent = 0

        # Add the thousands separator every 3 characters.
        if thousands_separators != '' and split > sep_range:
            # the first thousand separator may be located before 3 characters
            nmod = split % sep_range
            digs_b = digits[nmod:split]

            if nmod != 0:
                prev = digits[:nmod] + thousands_separators
            else:
                prev = ''

            dec_part = prev + thousands_separators.join(
                    digs_b[i:i+sep_range]
                    for i in range(0, split-nmod, sep_range)
                    )
        else:
            dec_part = digits[:split]

        # Finally, assemble the digits including the decimal point
        if precision == 0:
            return sign, dec_part + ('.' if alternate else '')

        digits = dec_part + "." + digits[split:]

    if strip_zeros:
        # Clean up trailing zeros
        digits = digits.rstrip('0')

    if digits[-1] == ".":
        if strip_last_zero:
            digits = digits[:-1]
        elif not alternate:
            digits += "0"

    return sign, digits


def format_scientific(s,
                      precision=6,
                      strip_zeros=False,
                      sign_spec='-',
                      base=10,
                      capitalize=False,
                      alternate=False):

    sep = 'E' if capitalize else 'e'

    # First, get the exponent to know how many digits we will need
    dps = precision+1
    sign, digits, exponent = to_digits_exp(
            s, max(dps+1, int(s[3]/blog2)), base)

    if sign != '-' and sign_spec != '-':
        sign = sign_spec

    # Rounding up kills some instances of "...99999"
    digits, exp_add = round_digits(digits, dps, base)
    exponent += exp_add

    if strip_zeros:
        # Clean up trailing zeros
        digits = digits.rstrip('0')
        precision = len(digits)

    if precision >= 1 and len(digits) > 1:
        return sign, digits[0] + '.' + digits[1:] + sep + f'{exponent:+03d}'
    else:
        if alternate:
            return sign, digits + '.' + sep + f'{exponent:+03d}'
        else:
            return sign, digits + sep + f'{exponent:+03d}'


def format_mpf(num, format_spec):
    format_dict = read_format_spec(format_spec)

    capitalize = False
    if format_dict['type'] in 'FGE':
        capitalize = True

    fmt_type = format_dict['type'].lower()
    precision = format_dict['precision']

    digits = ''
    sign = ''

    # Special cases:
    if num in (fzero, fnzero):
        if fmt_type in 'ef':
            digits = '0'
            if num[0]:
                sign = '-'

            if precision > 0:
                digits += '.' + precision*'0'
            else:
                if format_dict['alternate']:
                    digits += '.'

        elif fmt_type == 'g':
            if format_dict['alternate']:
                digits = '0.' + (precision-1)*'0'
            else:
                digits = '0'

        if fmt_type == 'e':
            sep = 'E' if capitalize else 'e'
            digits += sep + '+00'

    elif num == finf:
        if capitalize:
            digits = 'INF'
        else:
            digits = 'inf'

        if format_dict['sign'] == '+':
            digits = '+' + digits
    elif num == fninf:
        if capitalize:
            digits = '-INF'
        else:
            digits = '-inf'

    elif num == fnan:
        if capitalize:
            digits = 'NAN'
        else:
            digits = 'nan'

    # Now the general case
    else:
        strip_last_zero = False
        strip_zeros = False

        if fmt_type == 'g':
            if not format_dict['alternate']:
                strip_last_zero = True
                strip_zeros = True

            _, _, exp = to_digits_exp(num, 53/blog2, 10)
            if precision == 0:
                precision = 1

            if -4 <= exp < precision:
                fmt_type = 'f'
                precision = max(0, precision - exp - 1)
            else:
                fmt_type = 'e'
                precision = max(0, precision - 1)

        if fmt_type == 'f':
            sign, digits = format_fixed(
                    num,
                    precision=precision,
                    strip_zeros=strip_zeros,
                    strip_last_zero=strip_last_zero,
                    thousands_separators=format_dict['thousands_separators'],
                    sign_spec=format_dict['sign'],
                    sep_range=3,
                    base=10,
                    alternate=format_dict['alternate']
                    )
        elif fmt_type == 'e':
            sign, digits = format_scientific(
                    num,
                    precision=precision,
                    strip_zeros=strip_zeros,
                    sign_spec=format_dict['sign'],
                    base=10,
                    capitalize=capitalize,
                    alternate=format_dict['alternate']
                    )

        else:
            raise ValueError(
                    'Type {} not implemented yet'.format(format_dict['type']))

    nchars = len(digits) + len(sign)

    lpad, rpad = calc_padding(
            nchars, format_dict['width'], format_dict['align'])
    if format_dict['align'] == '=':
        return sign + lpad*format_dict['fill_char'] + digits + \
                rpad*format_dict['fill_char']
    else:
        return lpad*format_dict['fill_char'] + sign + digits \
                + rpad*format_dict['fill_char']
