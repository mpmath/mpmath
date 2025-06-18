import cmath
import ctypes
import math
import sys

import hypothesis.strategies as st
import pytest
from hypothesis import example, given, settings

from mpmath import fp, inf, mp, nan, ninf, workdps
from mpmath.libmp.libmpf import read_format_spec


vinfo = sys.version_info


@st.composite
def fmt_str(draw, types='fFeE', for_complex=False):
    res = ''

    # fill_char and align
    fill_char = draw(st.sampled_from(['']*3 + list('z;clxvjqwer')))
    if fill_char:
        skip_0_padding = True
        if for_complex:
            align = draw(st.sampled_from(list('<^>')))
        else:
            align = draw(st.sampled_from(list('<^>=')))
        res += fill_char + align
    else:
        align = draw(st.sampled_from([''] + list('<^>=')))
        if align == '=' and for_complex:
            align = ''
        if align:
            skip_0_padding = True
            res += align
        else:
            skip_0_padding = False

    # sign character
    res += draw(st.sampled_from([''] + list('-+ ')))

    # no_neg_0 (not used yet.)
    if vinfo >= (3, 11):
        res += draw(st.sampled_from([''] + ['z']))

    # alternate mode
    res += draw(st.sampled_from(['', '#']))

    # pad with 0s
    pad0 = draw(st.sampled_from(['', '0']))
    if pad0 and for_complex:
        pad0 = ''
    skip_thousand_separators = False
    if pad0 and not skip_0_padding:
        res += pad0
        skip_thousand_separators = True

    # Width
    res += draw(st.sampled_from(['']*7 + list(map(str, range(1, 40)))
                                + ([] if for_complex else ['0' + str(_)
                                                           for _ in range(40)])))

    # grouping character (thousand_separators)
    gchar = draw(st.sampled_from([''] + list(',_')))
    if gchar and not skip_thousand_separators:
        res += gchar

    # Precision
    prec = draw(st.sampled_from(['']*7 + list(map(str, range(40)))
                + ['0' + str(_) for _ in range(40)]))
    if prec:
        res += '.' + prec
        if vinfo >= (3, 14):
            gchar = draw(st.sampled_from([''] + list(',_')))
            res += gchar

    # Type
    res += draw(st.sampled_from(types))

    return res


def test_mpf_fmt_cpython():
    '''
    These tests assure that mpf.__format__ yields the same result as regular
    float.__format__, when dps is default.
    '''

    # 'f' code formatting.

    # zeros
    assert f'{mp.mpf(0):.0f}' == '0'
    assert f'{mp.mpf(0):.1f}' == '0.0'
    assert f'{mp.mpf(0):.2f}' == '0.00'
    assert f'{mp.mpf(0):.3f}' == '0.000'
    assert f'{mp.mpf(0):.50f}' == '0.00000000000000000000000000000000000000000000000000'

    # nan, infs
    assert f'{inf:f}' == 'inf'
    assert f'{inf:+f}' == '+inf'
    assert f'{inf:F}' == 'INF'
    assert f'{inf:+F}' == '+INF'

    assert f'{ninf:f}' == '-inf'
    assert f'{ninf:+f}' == '-inf'
    assert f'{ninf:F}' == '-INF'
    assert f'{ninf:+F}' == '-INF'

    assert f'{nan:f}' == 'nan'
    assert f'{nan:+f}' == '+nan'
    assert f'{nan:F}' == 'NAN'
    assert f'{nan:+F}' == '+NAN'

    # precision 0;  result should never include a .
    assert f'{mp.mpf(1.5):.0f}' == '2'
    assert f'{mp.mpf(2.5):.0f}' == '2'
    assert f'{mp.mpf(3.5):.0f}' == '4'
    assert f'{mp.mpf(0.0):.0f}' == '0'
    assert f'{mp.mpf(0.1):.0f}' == '0'
    assert f'{mp.mpf(0.001):.0f}' == '0'
    assert f'{mp.mpf(10.0):.0f}' == '10'
    assert f'{mp.mpf(10.1):.0f}' == '10'
    assert f'{mp.mpf(10.01):.0f}' == '10'
    assert f'{mp.mpf(123.456):.0f}' == '123'
    assert f'{mp.mpf(1234.56):.0f}' == '1235'
    assert f'{mp.mpf(1e49):.0f}' == '9999999999999999464902769475481793196872414789632'
    assert f'{mp.mpf(9.9999999999999987e+49):.0f}' == '99999999999999986860582406952576489172979654066176'
    assert f'{mp.mpf(1e50):.0f}' == '100000000000000007629769841091887003294964970946560'

    # precision 1
    assert f'{mp.mpf(0.0001):.1f}' == '0.0'
    assert f'{mp.mpf(0.001):.1f}' == '0.0'
    assert f'{mp.mpf(0.01):.1f}' == '0.0'
    assert f'{mp.mpf(0.04):.1f}' == '0.0'
    assert f'{mp.mpf(0.06):.1f}' == '0.1'
    assert f'{mp.mpf(0.25):.1f}' == '0.2'
    assert f'{mp.mpf(0.75):.1f}' == '0.8'
    assert f'{mp.mpf(1.4):.1f}' == '1.4'
    assert f'{mp.mpf(1.5):.1f}' == '1.5'
    assert f'{mp.mpf(10.0):.1f}' == '10.0'
    assert f'{mp.mpf(1000.03):.1f}' == '1000.0'
    assert f'{mp.mpf(1234.5678):.1f}' == '1234.6'
    assert f'{mp.mpf(1234.7499):.1f}' == '1234.7'
    assert f'{mp.mpf(1234.75):.1f}' == '1234.8'

    # precision 2
    assert f'{mp.mpf(0.0001):.2f}' == '0.00'
    assert f'{mp.mpf(0.001):.2f}' == '0.00'
    assert f'{mp.mpf(0.004999):.2f}' == '0.00'
    assert f'{mp.mpf(0.005001):.2f}' == '0.01'
    assert f'{mp.mpf(0.01):.2f}' == '0.01'
    assert f'{mp.mpf(0.125):.2f}' == '0.12'
    assert f'{mp.mpf(0.375):.2f}' == '0.38'
    assert f'{mp.mpf(1234500):.2f}' == '1234500.00'
    assert f'{mp.mpf(1234560):.2f}' == '1234560.00'
    assert f'{mp.mpf(1234567):.2f}' == '1234567.00'
    assert f'{mp.mpf(1234567.8):.2f}' == '1234567.80'
    assert f'{mp.mpf(1234567.89):.2f}' == '1234567.89'
    assert f'{mp.mpf(1234567.891):.2f}' == '1234567.89'
    assert f'{mp.mpf(1234567.8912):.2f}' == '1234567.89'

    # alternate form always includes a decimal point.  This only
    # makes a difference when the precision is 0.
    assert f'{mp.mpf(0):#.0f}' == '0.'
    assert f'{mp.mpf(0):#.1f}' == '0.0'
    assert f'{mp.mpf(1.5):#.0f}' == '2.'
    assert f'{mp.mpf(2.5):#.0f}' == '2.'
    assert f'{mp.mpf(10.1):#.0f}' == '10.'
    assert f'{mp.mpf(1234.56):#.0f}' == '1235.'
    assert f'{mp.mpf(1.4):#.1f}' == '1.4'
    assert f'{mp.mpf(0.375):#.2f}' == '0.38'

    # if precision is omitted it defaults to 6
    assert f'{mp.mpf(0):f}' == '0.000000'
    assert f'{mp.mpf(1230000):f}' == '1230000.000000'
    assert f'{mp.mpf(1234567):f}' == '1234567.000000'
    assert f'{mp.mpf(123.4567):f}' == '123.456700'
    assert f'{mp.mpf(1.23456789):f}' == '1.234568'
    assert f'{mp.mpf(0.00012):f}' == '0.000120'
    assert f'{mp.mpf(0.000123):f}' == '0.000123'
    assert f'{mp.mpf(0.00012345):f}' == '0.000123'
    assert f'{mp.mpf(0.000001):f}' == '0.000001'
    assert f'{mp.mpf(0.0000005001):f}' == '0.000001'
    assert f'{mp.mpf(0.0000004999):f}' == '0.000000'

    # grouping in fractional part
    assert f'{mp.mpf(0.0000004999):.9_f}' == '0.000_000_500'

    # 'e' code formatting with explicit precision (>= 0). Output should
    # always have exactly the number of places after the point that were
    # requested.

    # zeros
    assert f'{mp.mpf(0):.0e}' == '0e+00'
    assert f'{mp.mpf(0):.1e}' == '0.0e+00'
    assert f'{mp.mpf(0):.2e}' == '0.00e+00'
    assert f'{mp.mpf(0):.10e}' == '0.0000000000e+00'
    assert f'{mp.mpf(0):.50e}' == '0.00000000000000000000000000000000000000000000000000e+00'

    # nan, infs
    assert f'{inf:e}' == 'inf'
    assert f'{inf:+e}' == '+inf'
    assert f'{inf:E}' == 'INF'
    assert f'{inf:+E}' == '+INF'

    assert f'{ninf:e}' == '-inf'
    assert f'{ninf:+e}' == '-inf'
    assert f'{ninf:E}' == '-INF'
    assert f'{ninf:+E}' == '-INF'

    assert f'{nan:e}' == 'nan'
    assert f'{nan:+e}' == '+nan'
    assert f'{nan:E}' == 'NAN'
    assert f'{nan:+E}' == '+NAN'

    # precision 0.  no decimal point in the output
    assert f'{mp.mpf(0.01):.0e}' == '1e-02'
    assert f'{mp.mpf(0.1):.0e}' == '1e-01'
    assert f'{mp.mpf(1):.0e}' == '1e+00'
    assert f'{mp.mpf(10):.0e}' == '1e+01'
    assert f'{mp.mpf(100):.0e}' == '1e+02'
    assert f'{mp.mpf(0.012):.0e}' == '1e-02'
    assert f'{mp.mpf(0.12):.0e}' == '1e-01'
    assert f'{mp.mpf(1.2):.0e}' == '1e+00'
    assert f'{mp.mpf(12):.0e}' == '1e+01'
    assert f'{mp.mpf(120):.0e}' == '1e+02'
    assert f'{mp.mpf(123.456):.0e}' == '1e+02'
    assert f'{mp.mpf(0.000123456):.0e}' == '1e-04'
    assert f'{mp.mpf(123456000):.0e}' == '1e+08'
    assert f'{mp.mpf(0.5):.0e}' == '5e-01'
    assert f'{mp.mpf(1.4):.0e}' == '1e+00'
    assert f'{mp.mpf(1.5):.0e}' == '2e+00'
    assert f'{mp.mpf(1.6):.0e}' == '2e+00'
    assert f'{mp.mpf(2.4999999):.0e}' == '2e+00'
    assert f'{mp.mpf(2.5):.0e}' == '2e+00'
    assert f'{mp.mpf(2.5000001):.0e}' == '3e+00'
    assert f'{mp.mpf(3.499999999999):.0e}' == '3e+00'
    assert f'{mp.mpf(3.5):.0e}' == '4e+00'
    assert f'{mp.mpf(4.5):.0e}' == '4e+00'
    assert f'{mp.mpf(5.5):.0e}' == '6e+00'
    assert f'{mp.mpf(6.5):.0e}' == '6e+00'
    assert f'{mp.mpf(7.5):.0e}' == '8e+00'
    assert f'{mp.mpf(8.5):.0e}' == '8e+00'
    assert f'{mp.mpf(9.4999):.0e}' == '9e+00'
    assert f'{mp.mpf(9.5):.0e}' == '1e+01'
    assert f'{mp.mpf(10.5):.0e}' == '1e+01'
    assert f'{mp.mpf(14.999):.0e}' == '1e+01'
    assert f'{mp.mpf(15):.0e}' == '2e+01'

    # precision 1
    assert f'{mp.mpf(0.0001):.1e}' == '1.0e-04'
    assert f'{mp.mpf(0.001):.1e}' == '1.0e-03'
    assert f'{mp.mpf(0.01):.1e}' == '1.0e-02'
    assert f'{mp.mpf(0.1):.1e}' == '1.0e-01'
    assert f'{mp.mpf(1):.1e}' == '1.0e+00'
    assert f'{mp.mpf(10):.1e}' == '1.0e+01'
    assert f'{mp.mpf(100):.1e}' == '1.0e+02'
    assert f'{mp.mpf(120):.1e}' == '1.2e+02'
    assert f'{mp.mpf(123):.1e}' == '1.2e+02'
    assert f'{mp.mpf(123.4):.1e}' == '1.2e+02'

    # precision 2
    assert f'{mp.mpf(0.00013):.2e}' == '1.30e-04'
    assert f'{mp.mpf(0.000135):.2e}' == '1.35e-04'
    assert f'{mp.mpf(0.0001357):.2e}' == '1.36e-04'
    assert f'{mp.mpf(0.0001):.2e}' == '1.00e-04'
    assert f'{mp.mpf(0.001):.2e}' == '1.00e-03'
    assert f'{mp.mpf(0.01):.2e}' == '1.00e-02'
    assert f'{mp.mpf(0.1):.2e}' == '1.00e-01'
    assert f'{mp.mpf(1):.2e}' == '1.00e+00'
    assert f'{mp.mpf(10):.2e}' == '1.00e+01'
    assert f'{mp.mpf(100):.2e}' == '1.00e+02'
    assert f'{mp.mpf(1000):.2e}' == '1.00e+03'
    assert f'{mp.mpf(1500):.2e}' == '1.50e+03'
    assert f'{mp.mpf(1590):.2e}' == '1.59e+03'
    assert f'{mp.mpf(1598):.2e}' == '1.60e+03'
    assert f'{mp.mpf(1598.7):.2e}' == '1.60e+03'
    assert f'{mp.mpf(1598.76):.2e}' == '1.60e+03'
    assert f'{mp.mpf(9999):.2e}' == '1.00e+04'

    # omitted precision defaults to 6
    assert f'{mp.mpf(0):e}' == '0.000000e+00'
    assert f'{mp.mpf(165):e}' == '1.650000e+02'
    assert f'{mp.mpf(1234567):e}' == '1.234567e+06'
    assert f'{mp.mpf(12345678):e}' == '1.234568e+07'
    assert f'{mp.mpf(1.1):e}' == '1.100000e+00'

    # alternate form always contains a decimal point.  This only makes
    # a difference when precision is 0.

    assert f'{mp.mpf(0.01):#.0e}' == '1.e-02'
    assert f'{mp.mpf(0.1):#.0e}' == '1.e-01'
    assert f'{mp.mpf(1):#.0e}' == '1.e+00'
    assert f'{mp.mpf(10):#.0e}' == '1.e+01'
    assert f'{mp.mpf(100):#.0e}' == '1.e+02'
    assert f'{mp.mpf(0.012):#.0e}' == '1.e-02'
    assert f'{mp.mpf(0.12):#.0e}' == '1.e-01'
    assert f'{mp.mpf(1.2):#.0e}' == '1.e+00'
    assert f'{mp.mpf(12):#.0e}' == '1.e+01'
    assert f'{mp.mpf(120):#.0e}' == '1.e+02'
    assert f'{mp.mpf(123.456):#.0e}' == '1.e+02'
    assert f'{mp.mpf(0.000123456):#.0e}' == '1.e-04'
    assert f'{mp.mpf(123456000):#.0e}' == '1.e+08'
    assert f'{mp.mpf(0.5):#.0e}' == '5.e-01'
    assert f'{mp.mpf(1.4):#.0e}' == '1.e+00'
    assert f'{mp.mpf(1.5):#.0e}' == '2.e+00'
    assert f'{mp.mpf(1.6):#.0e}' == '2.e+00'
    assert f'{mp.mpf(2.4999999):#.0e}' == '2.e+00'
    assert f'{mp.mpf(2.5):#.0e}' == '2.e+00'
    assert f'{mp.mpf(2.5000001):#.0e}' == '3.e+00'
    assert f'{mp.mpf(3.499999999999):#.0e}' == '3.e+00'
    assert f'{mp.mpf(3.5):#.0e}' == '4.e+00'
    assert f'{mp.mpf(4.5):#.0e}' == '4.e+00'
    assert f'{mp.mpf(5.5):#.0e}' == '6.e+00'
    assert f'{mp.mpf(6.5):#.0e}' == '6.e+00'
    assert f'{mp.mpf(7.5):#.0e}' == '8.e+00'
    assert f'{mp.mpf(8.5):#.0e}' == '8.e+00'
    assert f'{mp.mpf(9.4999):#.0e}' == '9.e+00'
    assert f'{mp.mpf(9.5):#.0e}' == '1.e+01'
    assert f'{mp.mpf(10.5):#.0e}' == '1.e+01'
    assert f'{mp.mpf(14.999):#.0e}' == '1.e+01'
    assert f'{mp.mpf(15):#.0e}' == '2.e+01'
    assert f'{mp.mpf(123.4):#.1e}' == '1.2e+02'
    assert f'{mp.mpf(0.0001357):#.2e}' == '1.36e-04'

    # 'g' code formatting.

    # zeros
    assert f'{mp.mpf(0):.0g}' == '0'
    assert f'{mp.mpf(0):.1g}' == '0'
    assert f'{mp.mpf(0):.2g}' == '0'
    assert f'{mp.mpf(0):.3g}' == '0'
    assert f'{mp.mpf(0):.4g}' == '0'
    assert f'{mp.mpf(0):.10g}' == '0'
    assert f'{mp.mpf(0):.50g}' == '0'
    assert f'{mp.mpf(0):.100g}' == '0'

    # nan, infs
    assert f'{inf:g}' == 'inf'
    assert f'{inf:+g}' == '+inf'
    assert f'{inf:G}' == 'INF'
    assert f'{inf:+G}' == '+INF'

    assert f'{ninf:g}' == '-inf'
    assert f'{ninf:+g}' == '-inf'
    assert f'{ninf:G}' == '-INF'
    assert f'{ninf:+G}' == '-INF'

    assert f'{nan:g}' == 'nan'
    assert f'{nan:+g}' == '+nan'
    assert f'{nan:G}' == 'NAN'
    assert f'{nan:+G}' == '+NAN'

    # precision 0 doesn't make a lot of sense for the 'g' code (what does
    # it mean to have no significant digits?); in practice, it's interpreted
    # as identical to precision 1
    assert f'{mp.mpf(1000):.0g}' == '1e+03'
    assert f'{mp.mpf(100):.0g}' == '1e+02'
    assert f'{mp.mpf(10):.0g}' == '1e+01'
    assert f'{mp.mpf(1):.0g}' == '1'
    assert f'{mp.mpf(0.1):.0g}' == '0.1'
    assert f'{mp.mpf(0.01):.0g}' == '0.01'
    assert f'{mp.mpf(1e-3):.0g}' == '0.001'
    assert f'{mp.mpf(1e-4):.0g}' == '0.0001'
    assert f'{mp.mpf(1e-5):.0g}' == '1e-05'
    assert f'{mp.mpf(1e-6):.0g}' == '1e-06'
    assert f'{mp.mpf(12):.0g}' == '1e+01'
    assert f'{mp.mpf(120):.0g}' == '1e+02'
    assert f'{mp.mpf(1.2):.0g}' == '1'
    assert f'{mp.mpf(0.12):.0g}' == '0.1'
    assert f'{mp.mpf(0.012):.0g}' == '0.01'
    assert f'{mp.mpf(0.0012):.0g}' == '0.001'
    assert f'{mp.mpf(0.00012):.0g}' == '0.0001'
    assert f'{mp.mpf(0.000012):.0g}' == '1e-05'
    assert f'{mp.mpf(0.0000012):.0g}' == '1e-06'

    # precision 1 identical to precision 0
    assert f'{mp.mpf(1000):.1g}' == '1e+03'
    assert f'{mp.mpf(100):.1g}' == '1e+02'
    assert f'{mp.mpf(10):.1g}' == '1e+01'
    assert f'{mp.mpf(1):.1g}' == '1'
    assert f'{mp.mpf(0.1):.1g}' == '0.1'
    assert f'{mp.mpf(0.01):.1g}' == '0.01'
    assert f'{mp.mpf(1e-3):.1g}' == '0.001'
    assert f'{mp.mpf(1e-4):.1g}' == '0.0001'
    assert f'{mp.mpf(1e-5):.1g}' == '1e-05'
    assert f'{mp.mpf(1e-6):.1g}' == '1e-06'
    assert f'{mp.mpf(12):.1g}' == '1e+01'
    assert f'{mp.mpf(120):.1g}' == '1e+02'
    assert f'{mp.mpf(1.2):.1g}' == '1'
    assert f'{mp.mpf(0.12):.1g}' == '0.1'
    assert f'{mp.mpf(0.012):.1g}' == '0.01'
    assert f'{mp.mpf(0.0012):.1g}' == '0.001'
    assert f'{mp.mpf(0.00012):.1g}' == '0.0001'
    assert f'{mp.mpf(0.000012):.1g}' == '1e-05'
    assert f'{mp.mpf(0.0000012):.1g}' == '1e-06'

    # precision 2
    assert f'{mp.mpf(1000):.2g}' == '1e+03'
    assert f'{mp.mpf(100):.2g}' == '1e+02'
    assert f'{mp.mpf(10):.2g}' == '10'
    assert f'{mp.mpf(1):.2g}' == '1'
    assert f'{mp.mpf(0.1):.2g}' == '0.1'
    assert f'{mp.mpf(0.01):.2g}' == '0.01'
    assert f'{mp.mpf(0.001):.2g}' == '0.001'
    assert f'{mp.mpf(1e-4):.2g}' == '0.0001'
    assert f'{mp.mpf(1e-5):.2g}' == '1e-05'
    assert f'{mp.mpf(1e-6):.2g}' == '1e-06'
    assert f'{mp.mpf(1234):.2g}' == '1.2e+03'
    assert f'{mp.mpf(123):.2g}' == '1.2e+02'
    assert f'{mp.mpf(12.3):.2g}' == '12'
    assert f'{mp.mpf(1.23):.2g}' == '1.2'
    assert f'{mp.mpf(0.123):.2g}' == '0.12'
    assert f'{mp.mpf(0.0123):.2g}' == '0.012'
    assert f'{mp.mpf(0.00123):.2g}' == '0.0012'
    assert f'{mp.mpf(0.000123):.2g}' == '0.00012'
    assert f'{mp.mpf(0.0000123):.2g}' == '1.2e-05'

    # bad cases from http://bugs.python.org/issue9980
    assert f'{mp.mpf(38210.0):.12g}' == '38210'
    assert f'{mp.mpf(37210.0):.12g}' == '37210'
    assert f'{mp.mpf(36210.0):.12g}' == '36210'

    # alternate g formatting:  always include decimal point and
    # exactly <precision> significant digits.
    assert f'{mp.mpf(0):#.0g}' == '0.'
    assert f'{mp.mpf(0):#.1g}' == '0.'
    assert f'{mp.mpf(0):#.2g}' == '0.0'
    assert f'{mp.mpf(0):#.3g}' == '0.00'
    assert f'{mp.mpf(0):#.4g}' == '0.000'

    assert f'{mp.mpf(0.2):#.0g}' == '0.2'
    assert f'{mp.mpf(0.2):#.1g}' == '0.2'
    assert f'{mp.mpf(0.2):#.2g}' == '0.20'
    assert f'{mp.mpf(0.2):#.3g}' == '0.200'
    assert f'{mp.mpf(0.2):#.4g}' == '0.2000'
    assert f'{mp.mpf(0.2):#.10g}' == '0.2000000000'

    assert f'{mp.mpf(2):#.0g}' == '2.'
    assert f'{mp.mpf(2):#.1g}' == '2.'
    assert f'{mp.mpf(2):#.2g}' == '2.0'
    assert f'{mp.mpf(2):#.3g}' == '2.00'
    assert f'{mp.mpf(2):#.4g}' == '2.000'

    assert f'{mp.mpf(20):#.0g}' == '2.e+01'
    assert f'{mp.mpf(20):#.1g}' == '2.e+01'
    assert f'{mp.mpf(20):#.2g}' == '20.'
    assert f'{mp.mpf(20):#.3g}' == '20.0'
    assert f'{mp.mpf(20):#.4g}' == '20.00'

    assert f'{mp.mpf(234.56):#.0g}' == '2.e+02'
    assert f'{mp.mpf(234.56):#.1g}' == '2.e+02'
    assert f'{mp.mpf(234.56):#.2g}' == '2.3e+02'
    assert f'{mp.mpf(234.56):#.3g}' == '235.'
    assert f'{mp.mpf(234.56):#.4g}' == '234.6'
    assert f'{mp.mpf(234.56):#.5g}' == '234.56'
    assert f'{mp.mpf(234.56):#.6g}' == '234.560'

    # '%' code formatting.

    # nan, infs
    assert f'{inf:%}' == 'inf%'
    assert f'{inf:+%}' == '+inf%'

    assert f'{ninf:%}' == '-inf%'
    assert f'{ninf:+%}' == '-inf%'

    assert f'{nan:%}' == 'nan%'
    assert f'{nan:+%}' == '+nan%'

    # No formatting code.

    assert f'{mp.mpf(0.0):.0}' == '0e+00'
    mp.pretty_dps = 'repr'
    assert f'{mp.pi}' == '3.1415926535897931'
    mp.pretty_dps = 'str'
    assert f'{mp.pi}' == '3.14159265358979'


@settings(max_examples=20000)
@given(fmt_str(types=list('fFeEgG%') + ['']),
       st.floats(allow_nan=True,
                 allow_infinity=True,
                 allow_subnormal=True))
@example(fmt='.0g', x=9.995074823339339e-05)  # issue 880
@example(fmt='.016f', x=0.1)  # issue 915
@example(fmt='0030f', x=0.3)
@example(fmt='0=13,f', x=1.1)  # issue 917
@example(fmt='013,f', x=1.1)
@example(fmt='013,.0%', x=1.1)
@example(fmt='010.6,f', x=0.1234567891)
@example(fmt='010.7,f', x=0.1234567891)
@example(fmt='010._f', x=0.1234567891)
def test_mpf_floats_bulk(fmt, x):
    '''
    These are additional random tests that check that mp.mpf and fp.mpf yield
    the same results for default precision.
    '''

    mp.pretty_dps = "repr"
    if not x and math.copysign(1, x) == -1:
        return  # skip negative zero
    spec = read_format_spec(fmt)
    if spec['frac_separators'] and vinfo < (3, 14):
        mp.pretty_dps = "str"
        return  # see also python/cpython#130860
    if not spec['type'] and spec['precision'] < 0 and math.isfinite(x):
        # The mpmath could choose a different decimal
        # representative (wrt CPython) for same binary
        # floating-point number.
        assert float(format(x)) == float(format(mp.mpf(x)))
    else:
        if spec['type'] == '%' and math.isinf(100*x):
            return  # mpf can't overflow
        assert format(x, fmt) == format(mp.mpf(x), fmt)
    mp.pretty_dps = "str"


@settings(max_examples=20000)
@given(fmt_str(types=list('gGfFeE') + [''], for_complex=True),
       st.complex_numbers(allow_nan=True,
                          allow_infinity=True,
                          allow_subnormal=True))
def test_mpc_complexes(fmt, z):
    mp.pretty_dps = "repr"
    if ((not z.real and math.copysign(1, z.real) == -1)
            or (not z.imag and math.copysign(1, z.imag) == -1)):
        mp.pretty_dps = "str"
        return  # skip negative zero
    spec = read_format_spec(fmt)
    if spec['frac_separators'] and vinfo < (3, 14):
        return  # see also python/cpython#130860
    if spec['precision'] < 0 and any(math.isfinite(_) for _ in [z.real, z.imag]):
        # The mpmath could choose a different decimal
        # representative (wrt CPython) for same binary
        # floating-point number.
        if cmath.isnan(complex(format(z))):
            assert cmath.isnan(complex(format(mp.mpc(z))))
        else:
            assert complex(format(z)) == complex(format(mp.mpc(z)))
    else:
        assert format(z, fmt) == format(mp.mpc(z), fmt)
    mp.pretty_dps = "str"


def test_mpc_fmt():
    pytest.raises(ValueError, lambda: f'{mp.mpc(1j):=10f}')
    pytest.raises(ValueError, lambda: f'{mp.mpc(1j):010f}')
    pytest.raises(ValueError, lambda: f'{mp.mpc(1j):%}')

    assert f'{1+1j:.0g}' == f'{mp.mpc(1+1j):.0g}'
    assert f'{1+1.1j:.2g}' == f'{mp.mpc(1+1.1j):.2g}'


def test_mpf_fmt():
    '''
    These tests are either specific tests to mpf, or tests that cover
    code that is not covered in the CPython tests.
    '''

    with workdps(1000):
        # Numbers with more than 15 significant digits
        # fixed format
        assert f"{mp.mpf('1.234567890123456789'):.20f}" == '1.23456789012345678900'
        assert f"{mp.mpf('1.234567890123456789'):.25f}" == '1.2345678901234567890000000'
        assert f"{mp.mpf('1.234567890123456789'):.30f}" == '1.234567890123456789000000000000'
        assert f"{mp.mpf('1e-50'):.50f}" == '0.00000000000000000000000000000000000000000000000001'

        # scientific notation
        assert f"{mp.mpf('1.234567890123456789'):.20e}" == '1.23456789012345678900e+00'
        assert f"{mp.mpf('1.234567890123456789'):.25e}" == '1.2345678901234567890000000e+00'
        assert f"{mp.mpf('1.234567890123456789'):.30e}" == '1.234567890123456789000000000000e+00'
        assert f"{mp.mpf('1e-50'):.50e}" == '1.00000000000000000000000000000000000000000000000000e-50'

        # width and fill char
        assert f"{mp.mpf('0.01'):z<10.5f}" == '0.01000zzz'
        assert f"{mp.mpf('0.01'):z^10.5f}" == 'z0.01000zz'
        assert f"{mp.mpf('0.01'):z>10.5f}" == 'zzz0.01000'
        assert f"{mp.mpf('0.01'):z=10.5f}" == 'zzz0.01000'

        assert f"{mp.mpf('0.01'):z<+10.5f}" == '+0.01000zz'
        assert f"{mp.mpf('0.01'):z^+10.5f}" == 'z+0.01000z'
        assert f"{mp.mpf('0.01'):z>+10.5f}" == 'zz+0.01000'
        assert f"{mp.mpf('0.01'):z=+10.5f}" == '+zz0.01000'

        assert f"{mp.mpf('-0.01'):z<10.5f}" == '-0.01000zz'
        assert f"{mp.mpf('-0.01'):z^10.5f}" == 'z-0.01000z'
        assert f"{mp.mpf('-0.01'):z>10.5f}" == 'zz-0.01000'
        assert f"{mp.mpf('-0.01'):z=10.5f}" == '-zz0.01000'

        assert f"{mp.mpf('0.01'):z<15.5e}" == '1.00000e-02zzzz'
        assert f"{mp.mpf('0.01'):z^15.5e}" == 'zz1.00000e-02zz'
        assert f"{mp.mpf('0.01'):z>15.5e}" == 'zzzz1.00000e-02'
        assert f"{mp.mpf('0.01'):z=15.5e}" == 'zzzz1.00000e-02'

        assert f"{mp.mpf('0.01'):z<+15.5e}" == '+1.00000e-02zzz'
        assert f"{mp.mpf('0.01'):z^+15.5e}" == 'z+1.00000e-02zz'
        assert f"{mp.mpf('0.01'):z>+15.5e}" == 'zzz+1.00000e-02'
        assert f"{mp.mpf('0.01'):z=+15.5e}" == '+zzz1.00000e-02'

        assert f"{mp.mpf('-0.01'):z<15.5e}" == '-1.00000e-02zzz'
        assert f"{mp.mpf('-0.01'):z^15.5e}" == 'z-1.00000e-02zz'
        assert f"{mp.mpf('-0.01'):z>15.5e}" == 'zzz-1.00000e-02'
        assert f"{mp.mpf('-0.01'):z=15.5e}" == '-zzz1.00000e-02'

        # capitalized scientific notation
        assert f"{mp.mpf('-0.01'):z<15.5E}" == '-1.00000E-02zzz'

        # generalized format
        assert f"{mp.mpf('1.234567890123456789'):.20g}" == '1.234567890123456789'
        assert f"{mp.mpf('1.234567890123456789'):.25g}" == '1.234567890123456789'
        assert f"{mp.mpf('1.234567890123456789'):.30g}" == '1.234567890123456789'
        assert f"{mp.mpf('1e-50'):.50g}" == '1e-50'
        assert f"{mp.mpf('1e-50'):.50G}" == '1E-50'

        assert f"{mp.mpf('1e-51'):}" == '1e-51'

        # thousands separator
        assert f"{mp.mpf('1e9'):,.0f}" == '1,000,000,000'
        assert f"{mp.mpf('123456789.0123456'):,.4f}" == '123,456,789.0123'
        assert f"{mp.mpf('1234567890.123456'):_.4f}" == '1_234_567_890.1235'
        assert f"{mp.mpf('1234.5678'):_.4f}" == '1_234.5678'

        assert f"{mp.mpf('1e9'):,.0e}" == '1e+09'
        assert f"{mp.mpf('123456789.0123456'):,.4e}" == '1.2346e+08'
        assert f"{mp.mpf('1234567890.123456'):_.4e}" == '1.2346e+09'
        assert f"{mp.mpf('1234.5678'):_.4e}" == '1.2346e+03'

        # Tests for no_neg_0
        assert f"{mp.mpf('-1e-4'):,.2f}" == '-0.00'
        assert f"{mp.mpf('-1e-4'):z,.2f}" == '0.00'

        # Tests for = alignment
        assert f"{mp.mpf('0.24'):=+20.2f}" == '+               0.24'
        assert f"{mp.mpf('0.24'):=+020.2e}" == '+000000000002.40e-01'
        assert f"{mp.mpf('0.24'):=+020.2g}" == '+0000000000000000.24'

        # Tests for different kinds of rounding
        num = mp.mpf('-1.23456789999901234567')
        assert f"{num:=.2Uf}" == "-1.23"
        assert f"{num:=.2Df}" == "-1.24"
        assert f"{num:=.2Zf}" == "-1.23"
        assert f"{num:=.2Nf}" == "-1.23"
        assert f"{num:=.2Yf}" == "-1.24"

        assert f"{num:=.3Uf}" == "-1.234"
        assert f"{num:=.3Df}" == "-1.235"
        assert f"{num:=.3Zf}" == "-1.234"
        assert f"{num:=.3Nf}" == "-1.235"
        assert f"{num:=.3Yf}" == "-1.235"

        assert f"{num:=.10Uf}" == "-1.2345678999"
        assert f"{num:=.10Df}" == "-1.2345679000"
        assert f"{num:=.10Zf}" == "-1.2345678999"
        assert f"{num:=.10Nf}" == "-1.2345679000"
        assert f"{num:=.10Yf}" == "-1.2345679000"

        num = mp.mpf('1.23456789999901234567')
        assert f"{num:=.2Uf}" == "1.24"
        assert f"{num:=.2Df}" == "1.23"
        assert f"{num:=.2Zf}" == "1.23"
        assert f"{num:=.2Nf}" == "1.23"
        assert f"{num:=.2Yf}" == "1.24"

        assert f"{num:=.3Uf}" == "1.235"
        assert f"{num:=.3Df}" == "1.234"
        assert f"{num:=.3Zf}" == "1.234"
        assert f"{num:=.3Nf}" == "1.235"
        assert f"{num:=.3Yf}" == "1.235"

        assert f"{num:=.10Uf}" == "1.2345679000"
        assert f"{num:=.10Df}" == "1.2345678999"
        assert f"{num:=.10Zf}" == "1.2345678999"
        assert f"{num:=.10Nf}" == "1.2345679000"
        assert f"{num:=.10Yf}" == "1.2345679000"

        num = mp.mpf('-123.456789999901234567')
        assert f"{num:=.2Ue}" == "-1.23e+02"
        assert f"{num:=.2De}" == "-1.24e+02"
        assert f"{num:=.2Ze}" == "-1.23e+02"
        assert f"{num:=.2Ne}" == "-1.23e+02"
        assert f"{num:=.2Ye}" == "-1.24e+02"

        assert f"{num:=.3Ue}" == "-1.234e+02"
        assert f"{num:=.3De}" == "-1.235e+02"
        assert f"{num:=.3Ze}" == "-1.234e+02"
        assert f"{num:=.3Ne}" == "-1.235e+02"
        assert f"{num:=.3Ye}" == "-1.235e+02"

        assert f"{num:=.10Ue}" == "-1.2345678999e+02"
        assert f"{num:=.10De}" == "-1.2345679000e+02"
        assert f"{num:=.10Ze}" == "-1.2345678999e+02"
        assert f"{num:=.10Ne}" == "-1.2345679000e+02"
        assert f"{num:=.10Ye}" == "-1.2345679000e+02"

        num = mp.mpf('123456.789999901234567')
        assert f"{num:=.2Ue}" == "1.24e+05"
        assert f"{num:=.2De}" == "1.23e+05"
        assert f"{num:=.2Ze}" == "1.23e+05"
        assert f"{num:=.2Ne}" == "1.23e+05"
        assert f"{num:=.2Ye}" == "1.24e+05"

        assert f"{num:=.3Ue}" == "1.235e+05"
        assert f"{num:=.3De}" == "1.234e+05"
        assert f"{num:=.3Ze}" == "1.234e+05"
        assert f"{num:=.3Ne}" == "1.235e+05"
        assert f"{num:=.3Ye}" == "1.235e+05"

        assert f"{num:=.10Ue}" == "1.2345679000e+05"
        assert f"{num:=.10De}" == "1.2345678999e+05"
        assert f"{num:=.10Ze}" == "1.2345678999e+05"
        assert f"{num:=.10Ne}" == "1.2345679000e+05"
        assert f"{num:=.10Ye}" == "1.2345679000e+05"

        assert f"{mp.mpf('123.456'):.2Ug}" == "1.3e+02"
        assert f"{mp.mpf('123.456'):.2Dg}" == "1.2e+02"
        assert f"{mp.mpf('123.456'):.2Zg}" == "1.2e+02"
        assert f"{mp.mpf('123.456'):.2Ng}" == "1.2e+02"
        assert f"{mp.mpf('123.456'):.2Yg}" == "1.3e+02"

        assert f"{mp.mpf('-123.456'):.2Ug}" == "-1.2e+02"
        assert f"{mp.mpf('-123.456'):.2Dg}" == "-1.3e+02"
        assert f"{mp.mpf('-123.456'):.2Zg}" == "-1.2e+02"
        assert f"{mp.mpf('-123.456'):.2Ng}" == "-1.2e+02"
        assert f"{mp.mpf('-123.456'):.2Yg}" == "-1.3e+02"

        assert f"{mp.mpf('123.456'):.5Ug}" == "123.46"
        assert f"{mp.mpf('123.456'):.5Dg}" == "123.45"
        assert f"{mp.mpf('123.456'):.5Zg}" == "123.45"
        assert f"{mp.mpf('123.456'):.5Ng}" == "123.46"
        assert f"{mp.mpf('123.456'):.5Yg}" == "123.46"

        assert f"{mp.mpf('-123.456'):.5Ug}" == "-123.45"
        assert f"{mp.mpf('-123.456'):.5Dg}" == "-123.46"
        assert f"{mp.mpf('-123.456'):.5Zg}" == "-123.45"
        assert f"{mp.mpf('-123.456'):.5Ng}" == "-123.46"
        assert f"{mp.mpf('-123.456'):.5Yg}" == "-123.46"

        # Special cases were tying is relevant (cases involve exact floats)
        assert f"{mp.mpf('0.25'):.1Nf}" == "0.2"
        assert f"{mp.mpf('0.75'):.1Nf}" == "0.8"

        num = mp.mpf('0.1')
        assert f"{-num:=.2Df}" == "-0.11"
        assert f"{-num:=.3Df}" == "-0.101"
        assert f"{-num:=.4Df}" == "-0.1001"
        assert f"{-num:=.5Df}" == "-0.10001"
        assert f"{-num:=.6Df}" == "-0.100001"
        assert f"{-num:=.7Df}" == "-0.1000001"

        assert f"{-num:=.2De}" == "-1.01e-01"
        assert f"{-num:=.3De}" == "-1.001e-01"
        assert f"{-num:=.4De}" == "-1.0001e-01"
        assert f"{-num:=.5De}" == "-1.00001e-01"
        assert f"{-num:=.6De}" == "-1.000001e-01"
        assert f"{-num:=.7De}" == "-1.0000001e-01"

        assert f"{num:=.2Uf}" == "0.11"
        assert f"{num:=.3Uf}" == "0.101"
        assert f"{num:=.4Uf}" == "0.1001"
        assert f"{num:=.5Uf}" == "0.10001"
        assert f"{num:=.6Uf}" == "0.100001"
        assert f"{num:=.7Uf}" == "0.1000001"

        assert f"{num:=.2Ue}" == "1.01e-01"
        assert f"{num:=.3Ue}" == "1.001e-01"
        assert f"{num:=.4Ue}" == "1.0001e-01"
        assert f"{num:=.5Ue}" == "1.00001e-01"
        assert f"{num:=.6Ue}" == "1.000001e-01"
        assert f"{num:=.7Ue}" == "1.0000001e-01"

        num = mp.mpf('0.25')
        assert f"{-num:=.2Df}" == "-0.25"
        assert f"{-num:=.3Df}" == "-0.250"
        assert f"{-num:=.4Df}" == "-0.2500"
        assert f"{-num:=.5Df}" == "-0.25000"
        assert f"{-num:=.6Df}" == "-0.250000"
        assert f"{-num:=.7Df}" == "-0.2500000"

        assert f"{-num:=.2De}" == "-2.50e-01"
        assert f"{-num:=.3De}" == "-2.500e-01"
        assert f"{-num:=.4De}" == "-2.5000e-01"
        assert f"{-num:=.5De}" == "-2.50000e-01"
        assert f"{-num:=.6De}" == "-2.500000e-01"
        assert f"{-num:=.7De}" == "-2.5000000e-01"

        assert f"{num:=.2Uf}" == "0.25"
        assert f"{num:=.3Uf}" == "0.250"
        assert f"{num:=.4Uf}" == "0.2500"
        assert f"{num:=.5Uf}" == "0.25000"
        assert f"{num:=.6Uf}" == "0.250000"
        assert f"{num:=.7Uf}" == "0.2500000"

        assert f"{num:=.2Ue}" == "2.50e-01"
        assert f"{num:=.3Ue}" == "2.500e-01"
        assert f"{num:=.4Ue}" == "2.5000e-01"
        assert f"{num:=.5Ue}" == "2.50000e-01"
        assert f"{num:=.6Ue}" == "2.500000e-01"
        assert f"{num:=.7Ue}" == "2.5000000e-01"

        # Changing the work precision changes the number printed (This is
        # expected behavior)
        with mp.workdps(20):
            assert f"{mp.mpf('0.1'):=.4Uf}" == "0.1000"
            assert f"{mp.mpf('-0.1'):=.4Df}" == "-0.1000"


def test_default_rounding():
    x = mp.mpf(mp.pi)

    assert f"{x:.3f}" == '3.142'
    mp.rounding = 'd'
    assert f"{x:.3f}" == '3.141'
    mp.rounding = 'u'
    assert f"{x:.3f}" == '3.142'


def test_issue_858():
    for n in range(2, 15):
        str_num = '0.' + (n)*'9'
        fmt_str = '.' + str(n-1) + 'f'

        assert format(mp.mpf(str_num), fmt_str) == format(fp.mpf(str_num), fmt_str)

    assert format(mp.mpf(0.96875), '#.0g') == '1.'


def test_errors():
    with pytest.raises(ValueError):
        # wrong format type
        f"{mp.mpf('-4'):22.15k}"

    with pytest.raises(ValueError, match="Invalid format specifier '<z15.e'"):
        # no precision specified after .
        f"{mp.mpf('-0.01'):<z15.e}"

    with pytest.raises(ValueError, match="Invalid format specifier '10.5fk'"):
        f"{mp.mpf('4'):10.5fk}"

    with pytest.raises(ValueError, match="Invalid format specifier '12.3 E '"):
        f"{mp.mpf('4'):12.3 E }"

    with pytest.raises(ValueError):
        f"{mp.mpf(1):.f}"
    with pytest.raises(ValueError):
        f"{mp.mpf(1):._6f}"


@settings(max_examples=10000)
@given(st.floats(allow_nan=True, allow_infinity=True,
                 allow_subnormal=False))
@example(float('nan'))
@example(float('inf'))
def test_hexadecimal_bulk(x):
    if math.isnan(x):
        assert math.isnan(float.fromhex(f"{mp.mpf(x):a}"))
    else:
        assert float.fromhex(f"{mp.mpf(x):a}") == x


try:
    libc = ctypes.CDLL("libc.so.6")
except OSError:
    libc = None


def float_print(d, i):
    fmt = "%." + str(i) + "a\n"
    a = ctypes.create_string_buffer(256)
    libc.sprintf(a, bytes(fmt, 'utf-8'), ctypes.c_double(d))
    return a.raw.decode('utf-8').split("\n")[0]


@pytest.mark.skipif(libc is None, reason='requires libc')
@settings(max_examples=10000)
@given(st.floats(allow_nan=False, allow_infinity=False,
                 allow_subnormal=False),
       st.integers(min_value=0, max_value=15))
def test_hexadecimal_with_libc_bulk(x, p):
    fmt = '.' + str(p) + 'a'
    x_hex = float_print(x, p)
    m_hex = format(mp.mpf(x), fmt)
    assert mp.mpf(m_hex) == mp.mpf(x_hex)


@settings(max_examples=10000)
@given(st.floats(allow_nan=False, allow_infinity=False,
                 allow_subnormal=False),
       st.integers(min_value=-3, max_value=15))
def test_binary_with_gmpy2_bulk(x, p):
    gmpy2 = pytest.importorskip('gmpy2')
    if not x and math.copysign(1, x) == -1:
        return  # skip negative zero
    if p >= 0:
        fmt = '.' + str(p) + 'b'
    else:
        fmt = 'b'
    g_bin = format(gmpy2.mpfr(x), fmt)
    m_bin = format(mp.mpf(x), fmt)
    assert m_bin == g_bin


def test_binary_fmt():
    x = mp.mpf(3)
    assert f'{x:b}' == '1.1p+1'
    assert f'{x:.2b}' == '1.10p+1'
    assert f'{x:+.2b}' == '+1.10p+1'
    assert f'{x:#.2b}' == '0b1.10p+1'

    x = mp.mpf(0)
    assert f'{x:.2b}' == '0.00p+0'
    assert f'{x:b}' == '0p+0'


def test_hexadecimal_fmt():
    with workdps(1000):
        x = mp.mpf('1.234567890123456789')
        assert f'{x:.20a}' == '0x1.3c0ca428c59fb71a4194p+0'
        assert f'{x:.20A}' == '0X1.3C0CA428C59FB71A4194P+0'
        assert f'{x:.0a}' == '0x1p+0'
        assert f'{x:#.0a}' == '0x1.p+0'
        assert f"{mp.mpf('1.234567890123456789'):+.0a}" == '+0x1p+0'
