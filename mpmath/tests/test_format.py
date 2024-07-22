import random
import sys

import pytest

from mpmath import fp, inf, mp, nan, ninf, workdps


def random_fmt():
    '''
    This helper generates valid random format strings.
    '''

    fmt_str = '{:'

    # fill_char and align
    n = random.randint(0, 2)
    if n == 0:
        fmt_str += random.choice('z;clxvjqwer') + random.choice('<^>=')
        skip_0_padding = True
    elif n == 1:
        fmt_str += random.choice('<^>=')
        skip_0_padding = True
    else:
        skip_0_padding = False

    # sign character
    n = random.randint(0, 1)
    if n == 1:
        fmt_str += random.choice('-+ ')

    # no_neg_0 (not used yet.)
    if sys.version_info[:3] > (3, 11):
        n = random.randint(0, 1)
        if n == 1:
            fmt_str += 'z'

    # alternate mode
    n = random.randint(0, 1)
    if n == 1:
        fmt_str += '#'

    # pad with 0s
    n = random.randint(0, 1)
    skip_thousand_separators = False
    if n == 1 and not skip_0_padding:
        fmt_str += '0'
        skip_thousand_separators = True

    # Width
    n = random.randint(0, 2)
    if n > 0:
        fmt_str += str(random.randint(1, 40))

    # grouping character (thousand_separators)
    n = random.randint(0, 1)
    if n == 1 and not skip_thousand_separators:
        fmt_str += random.choice(',_')

    # Precision
    n = random.randint(0, 2)
    if n > 0:
        fmt_str += '.' + str(random.randint(1, 40))

    # Type
    fmt_str += random.choice('fFgGeE')
    fmt_str += '}'

    return fmt_str


def test_mpf_fmt_cpython():
    '''
    These tests assure that mpf.__format__ yields the same result as regular
    float.__format__, when dps is default.
    '''

    # zeros
    assert '{:.0f}'.format(mp.mpf(0)) == '0'
    assert '{:.1f}'.format(mp.mpf(0)) == '0.0'
    assert '{:.2f}'.format(mp.mpf(0)) == '0.00'
    assert '{:.3f}'.format(mp.mpf(0)) == '0.000'
    assert '{:.50f}'.format(mp.mpf(0)) == '0.00000000000000000000000000000000000000000000000000'

    # precision 0;  result should never include a .
    assert '{:.0f}'.format(mp.mpf(1.5)) == '2'
    assert '{:.0f}'.format(mp.mpf(2.5)) == '2'
    assert '{:.0f}'.format(mp.mpf(3.5)) == '4'
    assert '{:.0f}'.format(mp.mpf(0.0)) == '0'
    assert '{:.0f}'.format(mp.mpf(0.1)) == '0'
    assert '{:.0f}'.format(mp.mpf(0.001)) == '0'
    assert '{:.0f}'.format(mp.mpf(10.0)) == '10'
    assert '{:.0f}'.format(mp.mpf(10.1)) == '10'
    assert '{:.0f}'.format(mp.mpf(10.01)) == '10'
    assert '{:.0f}'.format(mp.mpf(123.456)) == '123'
    assert '{:.0f}'.format(mp.mpf(1234.56)) == '1235'
    assert '{:.0f}'.format(mp.mpf(1e49)) == '9999999999999999464902769475481793196872414789632'
    assert '{:.0f}'.format(mp.mpf(9.9999999999999987e+49)) == '99999999999999986860582406952576489172979654066176'
    assert '{:.0f}'.format(mp.mpf(1e50)) == '100000000000000007629769841091887003294964970946560'

    # precision 1
    assert '{:.1f}'.format(mp.mpf(0.0001)) == '0.0'
    assert '{:.1f}'.format(mp.mpf(0.001)) == '0.0'
    assert '{:.1f}'.format(mp.mpf(0.01)) == '0.0'
    assert '{:.1f}'.format(mp.mpf(0.04)) == '0.0'
    assert '{:.1f}'.format(mp.mpf(0.06)) == '0.1'
    assert '{:.1f}'.format(mp.mpf(0.25)) == '0.2'
    assert '{:.1f}'.format(mp.mpf(0.75)) == '0.8'
    assert '{:.1f}'.format(mp.mpf(1.4)) == '1.4'
    assert '{:.1f}'.format(mp.mpf(1.5)) == '1.5'
    assert '{:.1f}'.format(mp.mpf(10.0)) == '10.0'
    assert '{:.1f}'.format(mp.mpf(1000.03)) == '1000.0'
    assert '{:.1f}'.format(mp.mpf(1234.5678)) == '1234.6'
    assert '{:.1f}'.format(mp.mpf(1234.7499)) == '1234.7'
    assert '{:.1f}'.format(mp.mpf(1234.75)) == '1234.8'

    # precision 2
    assert '{:.2f}'.format(mp.mpf(0.0001)) == '0.00'
    assert '{:.2f}'.format(mp.mpf(0.001)) == '0.00'
    assert '{:.2f}'.format(mp.mpf(0.004999)) == '0.00'
    assert '{:.2f}'.format(mp.mpf(0.005001)) == '0.01'
    assert '{:.2f}'.format(mp.mpf(0.01)) == '0.01'
    assert '{:.2f}'.format(mp.mpf(0.125)) == '0.12'
    assert '{:.2f}'.format(mp.mpf(0.375)) == '0.38'
    assert '{:.2f}'.format(mp.mpf(1234500)) == '1234500.00'
    assert '{:.2f}'.format(mp.mpf(1234560)) == '1234560.00'
    assert '{:.2f}'.format(mp.mpf(1234567)) == '1234567.00'
    assert '{:.2f}'.format(mp.mpf(1234567.8)) == '1234567.80'
    assert '{:.2f}'.format(mp.mpf(1234567.89)) == '1234567.89'
    assert '{:.2f}'.format(mp.mpf(1234567.891)) == '1234567.89'
    assert '{:.2f}'.format(mp.mpf(1234567.8912)) == '1234567.89'

    # alternate form always includes a decimal point.  This only
    # makes a difference when the precision is 0.
    assert '{:#.0f}'.format(mp.mpf(0)) == '0.'
    assert '{:#.1f}'.format(mp.mpf(0)) == '0.0'
    assert '{:#.0f}'.format(mp.mpf(1.5)) == '2.'
    assert '{:#.0f}'.format(mp.mpf(2.5)) == '2.'
    assert '{:#.0f}'.format(mp.mpf(10.1)) == '10.'
    assert '{:#.0f}'.format(mp.mpf(1234.56)) == '1235.'
    assert '{:#.1f}'.format(mp.mpf(1.4)) == '1.4'
    assert '{:#.2f}'.format(mp.mpf(0.375)) == '0.38'

    # if precision is omitted it defaults to 6
    assert '{:f}'.format(mp.mpf(0)) == '0.000000'
    assert '{:f}'.format(mp.mpf(1230000)) == '1230000.000000'
    assert '{:f}'.format(mp.mpf(1234567)) == '1234567.000000'
    assert '{:f}'.format(mp.mpf(123.4567)) == '123.456700'
    assert '{:f}'.format(mp.mpf(1.23456789)) == '1.234568'
    assert '{:f}'.format(mp.mpf(0.00012)) == '0.000120'
    assert '{:f}'.format(mp.mpf(0.000123)) == '0.000123'
    assert '{:f}'.format(mp.mpf(0.00012345)) == '0.000123'
    assert '{:f}'.format(mp.mpf(0.000001)) == '0.000001'
    assert '{:f}'.format(mp.mpf(0.0000005001)) == '0.000001'
    assert '{:f}'.format(mp.mpf(0.0000004999)) == '0.000000'

    # 'e' code formatting with explicit precision (>= 0). Output should
    # always have exactly the number of places after the point that were
    # requested.

    # zeros
    assert '{:.0e}'.format(mp.mpf(0)) == '0e+00'
    assert '{:.1e}'.format(mp.mpf(0)) == '0.0e+00'
    assert '{:.2e}'.format(mp.mpf(0)) == '0.00e+00'
    assert '{:.10e}'.format(mp.mpf(0)) == '0.0000000000e+00'
    assert '{:.50e}'.format(mp.mpf(0)) == '0.00000000000000000000000000000000000000000000000000e+00'

    # precision 0.  no decimal point in the output
    assert '{:.0e}'.format(mp.mpf(0.01)) == '1e-02'
    assert '{:.0e}'.format(mp.mpf(0.1)) == '1e-01'
    assert '{:.0e}'.format(mp.mpf(1)) == '1e+00'
    assert '{:.0e}'.format(mp.mpf(10)) == '1e+01'
    assert '{:.0e}'.format(mp.mpf(100)) == '1e+02'
    assert '{:.0e}'.format(mp.mpf(0.012)) == '1e-02'
    assert '{:.0e}'.format(mp.mpf(0.12)) == '1e-01'
    assert '{:.0e}'.format(mp.mpf(1.2)) == '1e+00'
    assert '{:.0e}'.format(mp.mpf(12)) == '1e+01'
    assert '{:.0e}'.format(mp.mpf(120)) == '1e+02'
    assert '{:.0e}'.format(mp.mpf(123.456)) == '1e+02'
    assert '{:.0e}'.format(mp.mpf(0.000123456)) == '1e-04'
    assert '{:.0e}'.format(mp.mpf(123456000)) == '1e+08'
    assert '{:.0e}'.format(mp.mpf(0.5)) == '5e-01'
    assert '{:.0e}'.format(mp.mpf(1.4)) == '1e+00'
    assert '{:.0e}'.format(mp.mpf(1.5)) == '2e+00'
    assert '{:.0e}'.format(mp.mpf(1.6)) == '2e+00'
    assert '{:.0e}'.format(mp.mpf(2.4999999)) == '2e+00'
    assert '{:.0e}'.format(mp.mpf(2.5)) == '2e+00'
    assert '{:.0e}'.format(mp.mpf(2.5000001)) == '3e+00'
    assert '{:.0e}'.format(mp.mpf(3.499999999999)) == '3e+00'
    assert '{:.0e}'.format(mp.mpf(3.5)) == '4e+00'
    assert '{:.0e}'.format(mp.mpf(4.5)) == '4e+00'
    assert '{:.0e}'.format(mp.mpf(5.5)) == '6e+00'
    assert '{:.0e}'.format(mp.mpf(6.5)) == '6e+00'
    assert '{:.0e}'.format(mp.mpf(7.5)) == '8e+00'
    assert '{:.0e}'.format(mp.mpf(8.5)) == '8e+00'
    assert '{:.0e}'.format(mp.mpf(9.4999)) == '9e+00'
    assert '{:.0e}'.format(mp.mpf(9.5)) == '1e+01'
    assert '{:.0e}'.format(mp.mpf(10.5)) == '1e+01'
    assert '{:.0e}'.format(mp.mpf(14.999)) == '1e+01'
    assert '{:.0e}'.format(mp.mpf(15)) == '2e+01'

    # precision 1
    assert '{:.1e}'.format(mp.mpf(0.0001)) == '1.0e-04'
    assert '{:.1e}'.format(mp.mpf(0.001)) == '1.0e-03'
    assert '{:.1e}'.format(mp.mpf(0.01)) == '1.0e-02'
    assert '{:.1e}'.format(mp.mpf(0.1)) == '1.0e-01'
    assert '{:.1e}'.format(mp.mpf(1)) == '1.0e+00'
    assert '{:.1e}'.format(mp.mpf(10)) == '1.0e+01'
    assert '{:.1e}'.format(mp.mpf(100)) == '1.0e+02'
    assert '{:.1e}'.format(mp.mpf(120)) == '1.2e+02'
    assert '{:.1e}'.format(mp.mpf(123)) == '1.2e+02'
    assert '{:.1e}'.format(mp.mpf(123.4)) == '1.2e+02'

    # precision 2
    assert '{:.2e}'.format(mp.mpf(0.00013)) == '1.30e-04'
    assert '{:.2e}'.format(mp.mpf(0.000135)) == '1.35e-04'
    assert '{:.2e}'.format(mp.mpf(0.0001357)) == '1.36e-04'
    assert '{:.2e}'.format(mp.mpf(0.0001)) == '1.00e-04'
    assert '{:.2e}'.format(mp.mpf(0.001)) == '1.00e-03'
    assert '{:.2e}'.format(mp.mpf(0.01)) == '1.00e-02'
    assert '{:.2e}'.format(mp.mpf(0.1)) == '1.00e-01'
    assert '{:.2e}'.format(mp.mpf(1)) == '1.00e+00'
    assert '{:.2e}'.format(mp.mpf(10)) == '1.00e+01'
    assert '{:.2e}'.format(mp.mpf(100)) == '1.00e+02'
    assert '{:.2e}'.format(mp.mpf(1000)) == '1.00e+03'
    assert '{:.2e}'.format(mp.mpf(1500)) == '1.50e+03'
    assert '{:.2e}'.format(mp.mpf(1590)) == '1.59e+03'
    assert '{:.2e}'.format(mp.mpf(1598)) == '1.60e+03'
    assert '{:.2e}'.format(mp.mpf(1598.7)) == '1.60e+03'
    assert '{:.2e}'.format(mp.mpf(1598.76)) == '1.60e+03'
    assert '{:.2e}'.format(mp.mpf(9999)) == '1.00e+04'

    # omitted precision defaults to 6
    assert '{:e}'.format(mp.mpf(0)) == '0.000000e+00'
    assert '{:e}'.format(mp.mpf(165)) == '1.650000e+02'
    assert '{:e}'.format(mp.mpf(1234567)) == '1.234567e+06'
    assert '{:e}'.format(mp.mpf(12345678)) == '1.234568e+07'
    assert '{:e}'.format(mp.mpf(1.1)) == '1.100000e+00'

    # alternate form always contains a decimal point.  This only makes
    # a difference when precision is 0.

    assert '{:#.0e}'.format(mp.mpf(0.01)) == '1.e-02'
    assert '{:#.0e}'.format(mp.mpf(0.1)) == '1.e-01'
    assert '{:#.0e}'.format(mp.mpf(1)) == '1.e+00'
    assert '{:#.0e}'.format(mp.mpf(10)) == '1.e+01'
    assert '{:#.0e}'.format(mp.mpf(100)) == '1.e+02'
    assert '{:#.0e}'.format(mp.mpf(0.012)) == '1.e-02'
    assert '{:#.0e}'.format(mp.mpf(0.12)) == '1.e-01'
    assert '{:#.0e}'.format(mp.mpf(1.2)) == '1.e+00'
    assert '{:#.0e}'.format(mp.mpf(12)) == '1.e+01'
    assert '{:#.0e}'.format(mp.mpf(120)) == '1.e+02'
    assert '{:#.0e}'.format(mp.mpf(123.456)) == '1.e+02'
    assert '{:#.0e}'.format(mp.mpf(0.000123456)) == '1.e-04'
    assert '{:#.0e}'.format(mp.mpf(123456000)) == '1.e+08'
    assert '{:#.0e}'.format(mp.mpf(0.5)) == '5.e-01'
    assert '{:#.0e}'.format(mp.mpf(1.4)) == '1.e+00'
    assert '{:#.0e}'.format(mp.mpf(1.5)) == '2.e+00'
    assert '{:#.0e}'.format(mp.mpf(1.6)) == '2.e+00'
    assert '{:#.0e}'.format(mp.mpf(2.4999999)) == '2.e+00'
    assert '{:#.0e}'.format(mp.mpf(2.5)) == '2.e+00'
    assert '{:#.0e}'.format(mp.mpf(2.5000001)) == '3.e+00'
    assert '{:#.0e}'.format(mp.mpf(3.499999999999)) == '3.e+00'
    assert '{:#.0e}'.format(mp.mpf(3.5)) == '4.e+00'
    assert '{:#.0e}'.format(mp.mpf(4.5)) == '4.e+00'
    assert '{:#.0e}'.format(mp.mpf(5.5)) == '6.e+00'
    assert '{:#.0e}'.format(mp.mpf(6.5)) == '6.e+00'
    assert '{:#.0e}'.format(mp.mpf(7.5)) == '8.e+00'
    assert '{:#.0e}'.format(mp.mpf(8.5)) == '8.e+00'
    assert '{:#.0e}'.format(mp.mpf(9.4999)) == '9.e+00'
    assert '{:#.0e}'.format(mp.mpf(9.5)) == '1.e+01'
    assert '{:#.0e}'.format(mp.mpf(10.5)) == '1.e+01'
    assert '{:#.0e}'.format(mp.mpf(14.999)) == '1.e+01'
    assert '{:#.0e}'.format(mp.mpf(15)) == '2.e+01'
    assert '{:#.1e}'.format(mp.mpf(123.4)) == '1.2e+02'
    assert '{:#.2e}'.format(mp.mpf(0.0001357)) == '1.36e-04'

    # 'g' code formatting.

    # zeros
    assert '{:.0g}'.format(mp.mpf(0)) == '0'
    assert '{:.1g}'.format(mp.mpf(0)) == '0'
    assert '{:.2g}'.format(mp.mpf(0)) == '0'
    assert '{:.3g}'.format(mp.mpf(0)) == '0'
    assert '{:.4g}'.format(mp.mpf(0)) == '0'
    assert '{:.10g}'.format(mp.mpf(0)) == '0'
    assert '{:.50g}'.format(mp.mpf(0)) == '0'
    assert '{:.100g}'.format(mp.mpf(0)) == '0'

    # precision 0 doesn't make a lot of sense for the 'g' code (what does
    # it mean to have no significant digits?); in practice, it's interpreted
    # as identical to precision 1
    assert '{:.0g}'.format(mp.mpf(1000)) == '1e+03'
    assert '{:.0g}'.format(mp.mpf(100)) == '1e+02'
    assert '{:.0g}'.format(mp.mpf(10)) == '1e+01'
    assert '{:.0g}'.format(mp.mpf(1)) == '1'
    assert '{:.0g}'.format(mp.mpf(0.1)) == '0.1'
    assert '{:.0g}'.format(mp.mpf(0.01)) == '0.01'
    assert '{:.0g}'.format(mp.mpf(1e-3)) == '0.001'
    assert '{:.0g}'.format(mp.mpf(1e-4)) == '0.0001'
    assert '{:.0g}'.format(mp.mpf(1e-5)) == '1e-05'
    assert '{:.0g}'.format(mp.mpf(1e-6)) == '1e-06'
    assert '{:.0g}'.format(mp.mpf(12)) == '1e+01'
    assert '{:.0g}'.format(mp.mpf(120)) == '1e+02'
    assert '{:.0g}'.format(mp.mpf(1.2)) == '1'
    assert '{:.0g}'.format(mp.mpf(0.12)) == '0.1'
    assert '{:.0g}'.format(mp.mpf(0.012)) == '0.01'
    assert '{:.0g}'.format(mp.mpf(0.0012)) == '0.001'
    assert '{:.0g}'.format(mp.mpf(0.00012)) == '0.0001'
    assert '{:.0g}'.format(mp.mpf(0.000012)) == '1e-05'
    assert '{:.0g}'.format(mp.mpf(0.0000012)) == '1e-06'

    # precision 1 identical to precision 0
    assert '{:.1g}'.format(mp.mpf(1000)) == '1e+03'
    assert '{:.1g}'.format(mp.mpf(100)) == '1e+02'
    assert '{:.1g}'.format(mp.mpf(10)) == '1e+01'
    assert '{:.1g}'.format(mp.mpf(1)) == '1'
    assert '{:.1g}'.format(mp.mpf(0.1)) == '0.1'
    assert '{:.1g}'.format(mp.mpf(0.01)) == '0.01'
    assert '{:.1g}'.format(mp.mpf(1e-3)) == '0.001'
    assert '{:.1g}'.format(mp.mpf(1e-4)) == '0.0001'
    assert '{:.1g}'.format(mp.mpf(1e-5)) == '1e-05'
    assert '{:.1g}'.format(mp.mpf(1e-6)) == '1e-06'
    assert '{:.1g}'.format(mp.mpf(12)) == '1e+01'
    assert '{:.1g}'.format(mp.mpf(120)) == '1e+02'
    assert '{:.1g}'.format(mp.mpf(1.2)) == '1'
    assert '{:.1g}'.format(mp.mpf(0.12)) == '0.1'
    assert '{:.1g}'.format(mp.mpf(0.012)) == '0.01'
    assert '{:.1g}'.format(mp.mpf(0.0012)) == '0.001'
    assert '{:.1g}'.format(mp.mpf(0.00012)) == '0.0001'
    assert '{:.1g}'.format(mp.mpf(0.000012)) == '1e-05'
    assert '{:.1g}'.format(mp.mpf(0.0000012)) == '1e-06'

    # precision 2
    assert '{:.2g}'.format(mp.mpf(1000)) == '1e+03'
    assert '{:.2g}'.format(mp.mpf(100)) == '1e+02'
    assert '{:.2g}'.format(mp.mpf(10)) == '10'
    assert '{:.2g}'.format(mp.mpf(1)) == '1'
    assert '{:.2g}'.format(mp.mpf(0.1)) == '0.1'
    assert '{:.2g}'.format(mp.mpf(0.01)) == '0.01'
    assert '{:.2g}'.format(mp.mpf(0.001)) == '0.001'
    assert '{:.2g}'.format(mp.mpf(1e-4)) == '0.0001'
    assert '{:.2g}'.format(mp.mpf(1e-5)) == '1e-05'
    assert '{:.2g}'.format(mp.mpf(1e-6)) == '1e-06'
    assert '{:.2g}'.format(mp.mpf(1234)) == '1.2e+03'
    assert '{:.2g}'.format(mp.mpf(123)) == '1.2e+02'
    assert '{:.2g}'.format(mp.mpf(12.3)) == '12'
    assert '{:.2g}'.format(mp.mpf(1.23)) == '1.2'
    assert '{:.2g}'.format(mp.mpf(0.123)) == '0.12'
    assert '{:.2g}'.format(mp.mpf(0.0123)) == '0.012'
    assert '{:.2g}'.format(mp.mpf(0.00123)) == '0.0012'
    assert '{:.2g}'.format(mp.mpf(0.000123)) == '0.00012'
    assert '{:.2g}'.format(mp.mpf(0.0000123)) == '1.2e-05'

    # bad cases from http://bugs.python.org/issue9980
    assert '{:.12g}'.format(mp.mpf(38210.0)) == '38210'
    assert '{:.12g}'.format(mp.mpf(37210.0)) == '37210'
    assert '{:.12g}'.format(mp.mpf(36210.0)) == '36210'

    # alternate g formatting:  always include decimal point and
    # exactly <precision> significant digits.
    assert '{:#.0g}'.format(mp.mpf(0)) == '0.'
    assert '{:#.1g}'.format(mp.mpf(0)) == '0.'
    assert '{:#.2g}'.format(mp.mpf(0)) == '0.0'
    assert '{:#.3g}'.format(mp.mpf(0)) == '0.00'
    assert '{:#.4g}'.format(mp.mpf(0)) == '0.000'

    assert '{:#.0g}'.format(mp.mpf(0.2)) == '0.2'
    assert '{:#.1g}'.format(mp.mpf(0.2)) == '0.2'
    assert '{:#.2g}'.format(mp.mpf(0.2)) == '0.20'
    assert '{:#.3g}'.format(mp.mpf(0.2)) == '0.200'
    assert '{:#.4g}'.format(mp.mpf(0.2)) == '0.2000'
    assert '{:#.10g}'.format(mp.mpf(0.2)) == '0.2000000000'

    assert '{:#.0g}'.format(mp.mpf(2)) == '2.'
    assert '{:#.1g}'.format(mp.mpf(2)) == '2.'
    assert '{:#.2g}'.format(mp.mpf(2)) == '2.0'
    assert '{:#.3g}'.format(mp.mpf(2)) == '2.00'
    assert '{:#.4g}'.format(mp.mpf(2)) == '2.000'

    assert '{:#.0g}'.format(mp.mpf(20)) == '2.e+01'
    assert '{:#.1g}'.format(mp.mpf(20)) == '2.e+01'
    assert '{:#.2g}'.format(mp.mpf(20)) == '20.'
    assert '{:#.3g}'.format(mp.mpf(20)) == '20.0'
    assert '{:#.4g}'.format(mp.mpf(20)) == '20.00'

    assert '{:#.0g}'.format(mp.mpf(234.56)) == '2.e+02'
    assert '{:#.1g}'.format(mp.mpf(234.56)) == '2.e+02'
    assert '{:#.2g}'.format(mp.mpf(234.56)) == '2.3e+02'
    assert '{:#.3g}'.format(mp.mpf(234.56)) == '235.'
    assert '{:#.4g}'.format(mp.mpf(234.56)) == '234.6'
    assert '{:#.5g}'.format(mp.mpf(234.56)) == '234.56'
    assert '{:#.6g}'.format(mp.mpf(234.56)) == '234.560'


def test_mpf_float():
    '''
    These are additional random tests that check that mp.mpf and fp.mpf yield
    the same results for default precision.
    '''

    for _ in range(10000):
        fmt_str = random_fmt()
        num = random.uniform(-1e50, 1e50)

        assert fmt_str.format(fp.mpf(num)) == fmt_str.format(mp.mpf(num))


def test_mpf_fmt():
    '''
    These tests are specific to mpf.
    '''

    with workdps(1000):
        # Numbers with more than 15 significant digits
        # fixed format
        assert '{:.20f}'.format(mp.mpf('1.234567890123456789')) == '1.23456789012345678900'
        assert '{:.25f}'.format(mp.mpf('1.234567890123456789')) == '1.2345678901234567890000000'
        assert '{:.30f}'.format(mp.mpf('1.234567890123456789')) == '1.234567890123456789000000000000'
        assert '{:.50f}'.format(mp.mpf('1e-50')) == '0.00000000000000000000000000000000000000000000000001'

        # scientific notation
        assert '{:.20e}'.format(mp.mpf('1.234567890123456789')) == '1.23456789012345678900e+00'
        assert '{:.25e}'.format(mp.mpf('1.234567890123456789')) == '1.2345678901234567890000000e+00'
        assert '{:.30e}'.format(mp.mpf('1.234567890123456789')) == '1.234567890123456789000000000000e+00'
        assert '{:.50e}'.format(mp.mpf('1e-50')) == '1.00000000000000000000000000000000000000000000000000e-50'

        # width and fill char
        assert '{:z<10.5f}'.format(mp.mpf('0.01')) == '0.01000zzz'
        assert '{:z^10.5f}'.format(mp.mpf('0.01')) == 'z0.01000zz'
        assert '{:z>10.5f}'.format(mp.mpf('0.01')) == 'zzz0.01000'
        assert '{:z=10.5f}'.format(mp.mpf('0.01')) == 'zzz0.01000'

        assert '{:z<+10.5f}'.format(mp.mpf('0.01')) == '+0.01000zz'
        assert '{:z^+10.5f}'.format(mp.mpf('0.01')) == 'z+0.01000z'
        assert '{:z>+10.5f}'.format(mp.mpf('0.01')) == 'zz+0.01000'
        assert '{:z=+10.5f}'.format(mp.mpf('0.01')) == '+zz0.01000'

        assert '{:z<10.5f}'.format(mp.mpf('-0.01')) == '-0.01000zz'
        assert '{:z^10.5f}'.format(mp.mpf('-0.01')) == 'z-0.01000z'
        assert '{:z>10.5f}'.format(mp.mpf('-0.01')) == 'zz-0.01000'
        assert '{:z=10.5f}'.format(mp.mpf('-0.01')) == '-zz0.01000'

        assert '{:z<15.5e}'.format(mp.mpf('0.01')) == '1.00000e-02zzzz'
        assert '{:z^15.5e}'.format(mp.mpf('0.01')) == 'zz1.00000e-02zz'
        assert '{:z>15.5e}'.format(mp.mpf('0.01')) == 'zzzz1.00000e-02'
        assert '{:z=15.5e}'.format(mp.mpf('0.01')) == 'zzzz1.00000e-02'

        assert '{:z<+15.5e}'.format(mp.mpf('0.01')) == '+1.00000e-02zzz'
        assert '{:z^+15.5e}'.format(mp.mpf('0.01')) == 'z+1.00000e-02zz'
        assert '{:z>+15.5e}'.format(mp.mpf('0.01')) == 'zzz+1.00000e-02'
        assert '{:z=+15.5e}'.format(mp.mpf('0.01')) == '+zzz1.00000e-02'

        assert '{:z<15.5e}'.format(mp.mpf('-0.01')) == '-1.00000e-02zzz'
        assert '{:z^15.5e}'.format(mp.mpf('-0.01')) == 'z-1.00000e-02zz'
        assert '{:z>15.5e}'.format(mp.mpf('-0.01')) == 'zzz-1.00000e-02'
        assert '{:z=15.5e}'.format(mp.mpf('-0.01')) == '-zzz1.00000e-02'

        # capitalized scientific notation
        assert '{:z<15.5E}'.format(mp.mpf('-0.01')) == '-1.00000E-02zzz'

        # generalized format
        assert '{:.20g}'.format(mp.mpf('1.234567890123456789')) == '1.234567890123456789'
        assert '{:.25g}'.format(mp.mpf('1.234567890123456789')) == '1.234567890123456789'
        assert '{:.30g}'.format(mp.mpf('1.234567890123456789')) == '1.234567890123456789'
        assert '{:.50g}'.format(mp.mpf('1e-50')) == '1e-50'
        assert '{:.50G}'.format(mp.mpf('1e-50')) == '1E-50'

        # thousands separator
        assert '{:,.0f}'.format(mp.mpf('1e9')) == '1,000,000,000'
        assert '{:,.4f}'.format(mp.mpf('123456789.0123456')) == '123,456,789.0123'
        assert '{:_.4f}'.format(mp.mpf('1234567890.123456')) == '1_234_567_890.1235'
        assert '{:_.4f}'.format(mp.mpf('1234.5678')) == '1_234.5678'

        assert '{:,.0e}'.format(mp.mpf('1e9')) == '1e+09'
        assert '{:,.4e}'.format(mp.mpf('123456789.0123456')) == '1.2346e+08'
        assert '{:_.4e}'.format(mp.mpf('1234567890.123456')) == '1.2346e+09'
        assert '{:_.4e}'.format(mp.mpf('1234.5678')) == '1.2346e+03'


def test_mpf_fmt_special():
    assert '{:f}'.format(inf) == 'inf'
    assert '{:+f}'.format(inf) == '+inf'
    assert '{:F}'.format(inf) == 'INF'
    assert '{:+F}'.format(inf) == '+INF'

    assert '{:f}'.format(ninf) == '-inf'
    assert '{:+f}'.format(ninf) == '-inf'
    assert '{:F}'.format(ninf) == '-INF'
    assert '{:+F}'.format(ninf) == '-INF'

    assert '{:f}'.format(nan) == 'nan'
    assert '{:+f}'.format(nan) == '+nan'
    assert '{:F}'.format(nan) == 'NAN'
    assert '{:+F}'.format(nan) == '+NAN'

    assert '{:e}'.format(inf) == 'inf'
    assert '{:+e}'.format(inf) == '+inf'
    assert '{:E}'.format(inf) == 'INF'
    assert '{:+E}'.format(inf) == '+INF'

    assert '{:e}'.format(ninf) == '-inf'
    assert '{:+e}'.format(ninf) == '-inf'
    assert '{:E}'.format(ninf) == '-INF'
    assert '{:+E}'.format(ninf) == '-INF'

    assert '{:e}'.format(nan) == 'nan'
    assert '{:+e}'.format(nan) == '+nan'
    assert '{:E}'.format(nan) == 'NAN'
    assert '{:+E}'.format(nan) == '+NAN'

    assert '{:g}'.format(inf) == 'inf'
    assert '{:+g}'.format(inf) == '+inf'
    assert '{:G}'.format(inf) == 'INF'
    assert '{:+G}'.format(inf) == '+INF'

    assert '{:g}'.format(ninf) == '-inf'
    assert '{:+g}'.format(ninf) == '-inf'
    assert '{:G}'.format(ninf) == '-INF'
    assert '{:+G}'.format(ninf) == '-INF'

    assert '{:g}'.format(nan) == 'nan'
    assert '{:+g}'.format(nan) == '+nan'
    assert '{:G}'.format(nan) == 'NAN'
    assert '{:+G}'.format(nan) == '+NAN'

    # Test the same random formats with special numbers
    for num in (float('inf'), -float('inf'), float('nan'), 0):
        fmt_str = random_fmt()
        assert fmt_str.format(fp.mpf(num)) == fmt_str.format(mp.mpf(num))


def test_errors():
    with pytest.raises(Exception):
        # wrong format type
        assert '{:22.15k}'.format(mp.mpf('-4'))

    with pytest.raises(ValueError, match="Invalid format specifier '<z15.e'"):
        # no precision specified after .
        assert '{:<z15.e}'.format(mp.mpf('-0.01'))

    with pytest.raises(ValueError, match="Invalid format specifier '10.5fk'") as e:
        assert '{:10.5fk}'.format(mp.mpf('4'))

    with pytest.raises(ValueError, match="Invalid format specifier '12.3 E '") as e:
        assert '{:12.3 E }'.format(mp.mpf('4'))

    with pytest.raises(ValueError, match="Cannot specify both 0-padding "
                       "and a fill character") as e:
        assert '{:q<03f}'.format(mp.mpf('4'))
