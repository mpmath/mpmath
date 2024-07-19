import mpmath as mp

def test_mpf_fmt_cpython():
    '''
    These tests assure that mpf.__format__ yields the same result as regular
    float.__format__, when dps is default.
    '''
    # zeros
    assert '{:.0f}'.format(0) == '0'
    assert '{:.1f}'.format(0) == '0.0'
    assert '{:.2f}'.format(0) == '0.00'
    assert '{:.3f}'.format(0) == '0.000'
    assert '{:.50f}'.format(0) == '0.00000000000000000000000000000000000000000000000000'

    # precision 0;  result should never include a .
    assert '{:.0f}'.format(1.5) == '2'
    assert '{:.0f}'.format(2.5) == '2'
    assert '{:.0f}'.format(3.5) == '4'
    assert '{:.0f}'.format(0.0) == '0'
    assert '{:.0f}'.format(0.1) == '0'
    assert '{:.0f}'.format(0.001) == '0'
    assert '{:.0f}'.format(10.0) == '10'
    assert '{:.0f}'.format(10.1) == '10'
    assert '{:.0f}'.format(10.01) == '10'
    assert '{:.0f}'.format(123.456) == '123'
    assert '{:.0f}'.format(1234.56) == '1235'
    assert '{:.0f}'.format(1e49) == '9999999999999999464902769475481793196872414789632'
    assert '{:.0f}'.format(9.9999999999999987e+49) == '99999999999999986860582406952576489172979654066176'
    assert '{:.0f}'.format(1e50) == '100000000000000007629769841091887003294964970946560'

    # precision 1
    assert '{:.1f}'.format(0.0001) == '0.0'
    assert '{:.1f}'.format(0.001) == '0.0'
    assert '{:.1f}'.format(0.01) == '0.0'
    assert '{:.1f}'.format(0.04) == '0.0'
    assert '{:.1f}'.format(0.06) == '0.1'
    assert '{:.1f}'.format(0.25) == '0.2'
    assert '{:.1f}'.format(0.75) == '0.8'
    assert '{:.1f}'.format(1.4) == '1.4'
    assert '{:.1f}'.format(1.5) == '1.5'
    assert '{:.1f}'.format(10.0) == '10.0'
    assert '{:.1f}'.format(1000.03) == '1000.0'
    assert '{:.1f}'.format(1234.5678) == '1234.6'
    assert '{:.1f}'.format(1234.7499) == '1234.7'
    assert '{:.1f}'.format(1234.75) == '1234.8'

    # precision 2
    assert '{:.2f}'.format(0.0001) == '0.00'
    assert '{:.2f}'.format(0.001) == '0.00'
    assert '{:.2f}'.format(0.004999) == '0.00'
    assert '{:.2f}'.format(0.005001) == '0.01'
    assert '{:.2f}'.format(0.01) == '0.01'
    assert '{:.2f}'.format(0.125) == '0.12'
    assert '{:.2f}'.format(0.375) == '0.38'
    assert '{:.2f}'.format(1234500) == '1234500.00'
    assert '{:.2f}'.format(1234560) == '1234560.00'
    assert '{:.2f}'.format(1234567) == '1234567.00'
    assert '{:.2f}'.format(1234567.8) == '1234567.80'
    assert '{:.2f}'.format(1234567.89) == '1234567.89'
    assert '{:.2f}'.format(1234567.891) == '1234567.89'
    assert '{:.2f}'.format(1234567.8912) == '1234567.89'

    # alternate form always includes a decimal point.  This only
    # makes a difference when the precision is 0.
    assert '{:#.0f}'.format(0) == '0.'
    assert '{:#.1f}'.format(0) == '0.0'
    assert '{:#.0f}'.format(1.5) == '2.'
    assert '{:#.0f}'.format(2.5) == '2.'
    assert '{:#.0f}'.format(10.1) == '10.'
    assert '{:#.0f}'.format(1234.56) == '1235.'
    assert '{:#.1f}'.format(1.4) == '1.4'
    assert '{:#.2f}'.format(0.375) == '0.38'

    # if precision is omitted it defaults to 6
    assert '{:f}'.format(0) == '0.000000'
    assert '{:f}'.format(1230000) == '1230000.000000'
    assert '{:f}'.format(1234567) == '1234567.000000'
    assert '{:f}'.format(123.4567) == '123.456700'
    assert '{:f}'.format(1.23456789) == '1.234568'
    assert '{:f}'.format(0.00012) == '0.000120'
    assert '{:f}'.format(0.000123) == '0.000123'
    assert '{:f}'.format(0.00012345) == '0.000123'
    assert '{:f}'.format(0.000001) == '0.000001'
    assert '{:f}'.format(0.0000005001) == '0.000001'
    assert '{:f}'.format(0.0000004999) == '0.000000'

    # 'e' code formatting with explicit precision (>= 0). Output should
    # always have exactly the number of places after the point that were
    # requested.

    # zeros
    assert '{:.0e}'.format(0) == '0e+00'
    assert '{:.1e}'.format(0) == '0.0e+00'
    assert '{:.2e}'.format(0) == '0.00e+00'
    assert '{:.10e}'.format(0) == '0.0000000000e+00'
    assert '{:.50e}'.format(0) == '0.00000000000000000000000000000000000000000000000000e+00'

    # precision 0.  no decimal point in the output
    assert '{:.0e}'.format(0.01) == '1e-02'
    assert '{:.0e}'.format(0.1) == '1e-01'
    assert '{:.0e}'.format(1) == '1e+00'
    assert '{:.0e}'.format(10) == '1e+01'
    assert '{:.0e}'.format(100) == '1e+02'
    assert '{:.0e}'.format(0.012) == '1e-02'
    assert '{:.0e}'.format(0.12) == '1e-01'
    assert '{:.0e}'.format(1.2) == '1e+00'
    assert '{:.0e}'.format(12) == '1e+01'
    assert '{:.0e}'.format(120) == '1e+02'
    assert '{:.0e}'.format(123.456) == '1e+02'
    assert '{:.0e}'.format(0.000123456) == '1e-04'
    assert '{:.0e}'.format(123456000) == '1e+08'
    assert '{:.0e}'.format(0.5) == '5e-01'
    assert '{:.0e}'.format(1.4) == '1e+00'
    assert '{:.0e}'.format(1.5) == '2e+00'
    assert '{:.0e}'.format(1.6) == '2e+00'
    assert '{:.0e}'.format(2.4999999) == '2e+00'
    assert '{:.0e}'.format(2.5) == '2e+00'
    assert '{:.0e}'.format(2.5000001) == '3e+00'
    assert '{:.0e}'.format(3.499999999999) == '3e+00'
    assert '{:.0e}'.format(3.5) == '4e+00'
    assert '{:.0e}'.format(4.5) == '4e+00'
    assert '{:.0e}'.format(5.5) == '6e+00'
    assert '{:.0e}'.format(6.5) == '6e+00'
    assert '{:.0e}'.format(7.5) == '8e+00'
    assert '{:.0e}'.format(8.5) == '8e+00'
    assert '{:.0e}'.format(9.4999) == '9e+00'
    assert '{:.0e}'.format(9.5) == '1e+01'
    assert '{:.0e}'.format(10.5) == '1e+01'
    assert '{:.0e}'.format(14.999) == '1e+01'
    assert '{:.0e}'.format(15) == '2e+01'

    # precision 1
    assert '{:.1e}'.format(0.0001) == '1.0e-04'
    assert '{:.1e}'.format(0.001) == '1.0e-03'
    assert '{:.1e}'.format(0.01) == '1.0e-02'
    assert '{:.1e}'.format(0.1) == '1.0e-01'
    assert '{:.1e}'.format(1) == '1.0e+00'
    assert '{:.1e}'.format(10) == '1.0e+01'
    assert '{:.1e}'.format(100) == '1.0e+02'
    assert '{:.1e}'.format(120) == '1.2e+02'
    assert '{:.1e}'.format(123) == '1.2e+02'
    assert '{:.1e}'.format(123.4) == '1.2e+02'

    # precision 2
    assert '{:.2e}'.format(0.00013) == '1.30e-04'
    assert '{:.2e}'.format(0.000135) == '1.35e-04'
    assert '{:.2e}'.format(0.0001357) == '1.36e-04'
    assert '{:.2e}'.format(0.0001) == '1.00e-04'
    assert '{:.2e}'.format(0.001) == '1.00e-03'
    assert '{:.2e}'.format(0.01) == '1.00e-02'
    assert '{:.2e}'.format(0.1) == '1.00e-01'
    assert '{:.2e}'.format(1) == '1.00e+00'
    assert '{:.2e}'.format(10) == '1.00e+01'
    assert '{:.2e}'.format(100) == '1.00e+02'
    assert '{:.2e}'.format(1000) == '1.00e+03'
    assert '{:.2e}'.format(1500) == '1.50e+03'
    assert '{:.2e}'.format(1590) == '1.59e+03'
    assert '{:.2e}'.format(1598) == '1.60e+03'
    assert '{:.2e}'.format(1598.7) == '1.60e+03'
    assert '{:.2e}'.format(1598.76) == '1.60e+03'
    assert '{:.2e}'.format(9999) == '1.00e+04'

    # omitted precision defaults to 6
    assert '{:e}'.format(0) == '0.000000e+00'
    assert '{:e}'.format(165) == '1.650000e+02'
    assert '{:e}'.format(1234567) == '1.234567e+06'
    assert '{:e}'.format(12345678) == '1.234568e+07'
    assert '{:e}'.format(1.1) == '1.100000e+00'

    # alternate form always contains a decimal point.  This only makes
    # a difference when precision is 0.

    assert '{:#.0e}'.format(0.01) == '1.e-02'
    assert '{:#.0e}'.format(0.1) == '1.e-01'
    assert '{:#.0e}'.format(1) == '1.e+00'
    assert '{:#.0e}'.format(10) == '1.e+01'
    assert '{:#.0e}'.format(100) == '1.e+02'
    assert '{:#.0e}'.format(0.012) == '1.e-02'
    assert '{:#.0e}'.format(0.12) == '1.e-01'
    assert '{:#.0e}'.format(1.2) == '1.e+00'
    assert '{:#.0e}'.format(12) == '1.e+01'
    assert '{:#.0e}'.format(120) == '1.e+02'
    assert '{:#.0e}'.format(123.456) == '1.e+02'
    assert '{:#.0e}'.format(0.000123456) == '1.e-04'
    assert '{:#.0e}'.format(123456000) == '1.e+08'
    assert '{:#.0e}'.format(0.5) == '5.e-01'
    assert '{:#.0e}'.format(1.4) == '1.e+00'
    assert '{:#.0e}'.format(1.5) == '2.e+00'
    assert '{:#.0e}'.format(1.6) == '2.e+00'
    assert '{:#.0e}'.format(2.4999999) == '2.e+00'
    assert '{:#.0e}'.format(2.5) == '2.e+00'
    assert '{:#.0e}'.format(2.5000001) == '3.e+00'
    assert '{:#.0e}'.format(3.499999999999) == '3.e+00'
    assert '{:#.0e}'.format(3.5) == '4.e+00'
    assert '{:#.0e}'.format(4.5) == '4.e+00'
    assert '{:#.0e}'.format(5.5) == '6.e+00'
    assert '{:#.0e}'.format(6.5) == '6.e+00'
    assert '{:#.0e}'.format(7.5) == '8.e+00'
    assert '{:#.0e}'.format(8.5) == '8.e+00'
    assert '{:#.0e}'.format(9.4999) == '9.e+00'
    assert '{:#.0e}'.format(9.5) == '1.e+01'
    assert '{:#.0e}'.format(10.5) == '1.e+01'
    assert '{:#.0e}'.format(14.999) == '1.e+01'
    assert '{:#.0e}'.format(15) == '2.e+01'
    assert '{:#.1e}'.format(123.4) == '1.2e+02'
    assert '{:#.2e}'.format(0.0001357) == '1.36e-04'

    # 'g' code formatting.

    # zeros
    assert '{:.0g}'.format(0) == '0'
    assert '{:.1g}'.format(0) == '0'
    assert '{:.2g}'.format(0) == '0'
    assert '{:.3g}'.format(0) == '0'
    assert '{:.4g}'.format(0) == '0'
    assert '{:.10g}'.format(0) == '0'
    assert '{:.50g}'.format(0) == '0'
    assert '{:.100g}'.format(0) == '0'

    # precision 0 doesn't make a lot of sense for the 'g' code (what does
    # it mean to have no significant digits?); in practice, it's interpreted
    # as identical to precision 1
    assert '{:.0g}'.format(1000) == '1e+03'
    assert '{:.0g}'.format(100) == '1e+02'
    assert '{:.0g}'.format(10) == '1e+01'
    assert '{:.0g}'.format(1) == '1'
    assert '{:.0g}'.format(0.1) == '0.1'
    assert '{:.0g}'.format(0.01) == '0.01'
    assert '{:.0g}'.format(1e-3) == '0.001'
    assert '{:.0g}'.format(1e-4) == '0.0001'
    assert '{:.0g}'.format(1e-5) == '1e-05'
    assert '{:.0g}'.format(1e-6) == '1e-06'
    assert '{:.0g}'.format(12) == '1e+01'
    assert '{:.0g}'.format(120) == '1e+02'
    assert '{:.0g}'.format(1.2) == '1'
    assert '{:.0g}'.format(0.12) == '0.1'
    assert '{:.0g}'.format(0.012) == '0.01'
    assert '{:.0g}'.format(0.0012) == '0.001'
    assert '{:.0g}'.format(0.00012) == '0.0001'
    assert '{:.0g}'.format(0.000012) == '1e-05'
    assert '{:.0g}'.format(0.0000012) == '1e-06'

    # precision 1 identical to precision 0
    assert '{:.1g}'.format(1000) == '1e+03'
    assert '{:.1g}'.format(100) == '1e+02'
    assert '{:.1g}'.format(10) == '1e+01'
    assert '{:.1g}'.format(1) == '1'
    assert '{:.1g}'.format(0.1) == '0.1'
    assert '{:.1g}'.format(0.01) == '0.01'
    assert '{:.1g}'.format(1e-3) == '0.001'
    assert '{:.1g}'.format(1e-4) == '0.0001'
    assert '{:.1g}'.format(1e-5) == '1e-05'
    assert '{:.1g}'.format(1e-6) == '1e-06'
    assert '{:.1g}'.format(12) == '1e+01'
    assert '{:.1g}'.format(120) == '1e+02'
    assert '{:.1g}'.format(1.2) == '1'
    assert '{:.1g}'.format(0.12) == '0.1'
    assert '{:.1g}'.format(0.012) == '0.01'
    assert '{:.1g}'.format(0.0012) == '0.001'
    assert '{:.1g}'.format(0.00012) == '0.0001'
    assert '{:.1g}'.format(0.000012) == '1e-05'
    assert '{:.1g}'.format(0.0000012) == '1e-06'

    # precision 2
    assert '{:.2g}'.format(1000) == '1e+03'
    assert '{:.2g}'.format(100) == '1e+02'
    assert '{:.2g}'.format(10) == '10'
    assert '{:.2g}'.format(1) == '1'
    assert '{:.2g}'.format(0.1) == '0.1'
    assert '{:.2g}'.format(0.01) == '0.01'
    assert '{:.2g}'.format(0.001) == '0.001'
    assert '{:.2g}'.format(1e-4) == '0.0001'
    assert '{:.2g}'.format(1e-5) == '1e-05'
    assert '{:.2g}'.format(1e-6) == '1e-06'
    assert '{:.2g}'.format(1234) == '1.2e+03'
    assert '{:.2g}'.format(123) == '1.2e+02'
    assert '{:.2g}'.format(12.3) == '12'
    assert '{:.2g}'.format(1.23) == '1.2'
    assert '{:.2g}'.format(0.123) == '0.12'
    assert '{:.2g}'.format(0.0123) == '0.012'
    assert '{:.2g}'.format(0.00123) == '0.0012'
    assert '{:.2g}'.format(0.000123) == '0.00012'
    assert '{:.2g}'.format(0.0000123) == '1.2e-05'

    # bad cases from http://bugs.python.org/issue9980
    assert '{:.12g}'.format(38210.0) == '38210'
    assert '{:.12g}'.format(37210.0) == '37210'
    assert '{:.12g}'.format(36210.0) == '36210'

    # alternate g formatting:  always include decimal point and
    # exactly <precision> significant digits.
    assert '{:#.0g}'.format(0) == '0.'
    assert '{:#.1g}'.format(0) == '0.'
    assert '{:#.2g}'.format(0) == '0.0'
    assert '{:#.3g}'.format(0) == '0.00'
    assert '{:#.4g}'.format(0) == '0.000'

    assert '{:#.0g}'.format(0.2) == '0.2'
    assert '{:#.1g}'.format(0.2) == '0.2'
    assert '{:#.2g}'.format(0.2) == '0.20'
    assert '{:#.3g}'.format(0.2) == '0.200'
    assert '{:#.4g}'.format(0.2) == '0.2000'
    assert '{:#.10g}'.format(0.2) == '0.2000000000'

    assert '{:#.0g}'.format(2) == '2.'
    assert '{:#.1g}'.format(2) == '2.'
    assert '{:#.2g}'.format(2) == '2.0'
    assert '{:#.3g}'.format(2) == '2.00'
    assert '{:#.4g}'.format(2) == '2.000'

    assert '{:#.0g}'.format(20) == '2.e+01'
    assert '{:#.1g}'.format(20) == '2.e+01'
    assert '{:#.2g}'.format(20) == '20.'
    assert '{:#.3g}'.format(20) == '20.0'
    assert '{:#.4g}'.format(20) == '20.00'

    assert '{:#.0g}'.format(234.56) == '2.e+02'
    assert '{:#.1g}'.format(234.56) == '2.e+02'
    assert '{:#.2g}'.format(234.56) == '2.3e+02'
    assert '{:#.3g}'.format(234.56) == '235.'
    assert '{:#.4g}'.format(234.56) == '234.6'
    assert '{:#.5g}'.format(234.56) == '234.56'
    assert '{:#.6g}'.format(234.56) == '234.560'
    #with open(format_testfile, encoding="utf-8") as testfile:
    #    for line in testfile:
    #        if line.startswith('--') or line.startswith('%r'):
    #            continue
    #        line = line.strip()
    #        if not line:
    #            continue

    #        lhs, rhs = map(str.strip, line.split('->'))
    #        fmt, arg = lhs.split()
    #        f = mp.mpf(arg)

    #        if fmt != '%r':
    #            fmt2 = fmt[1:]
    #            assert format(f, fmt2) == rhs
    #            # Negative 0 is not covered by mpmath

    #            if float(rhs) != 0:
    #                assert format(-f, fmt2) == '-' + rhs


def test_mpf_fmt():
    '''
    These tests are specific to mpf.
    '''

    with mp.workdps(1000):
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
    inf = mp.inf
    ninf = mp.ninf
    nan = mp.nan

    assert '{:f}'.format(inf) == 'inf'
    assert '{:+f}'.format(inf) == '+inf'
    assert '{:F}'.format(inf) == 'INF'
    assert '{:+F}'.format(inf) == '+INF'

    assert '{:f}'.format(ninf) == '-inf'
    assert '{:+f}'.format(ninf) == '-inf'
    assert '{:F}'.format(ninf) == '-INF'
    assert '{:+F}'.format(ninf) == '-INF'

    assert '{:f}'.format(nan) == 'nan'
    assert '{:+f}'.format(nan) == 'nan'
    assert '{:F}'.format(nan) == 'NAN'
    assert '{:+F}'.format(nan) == 'NAN'

    assert '{:e}'.format(inf) == 'inf'
    assert '{:+e}'.format(inf) == '+inf'
    assert '{:E}'.format(inf) == 'INF'
    assert '{:+E}'.format(inf) == '+INF'

    assert '{:e}'.format(ninf) == '-inf'
    assert '{:+e}'.format(ninf) == '-inf'
    assert '{:E}'.format(ninf) == '-INF'
    assert '{:+E}'.format(ninf) == '-INF'

    assert '{:e}'.format(nan) == 'nan'
    assert '{:+e}'.format(nan) == 'nan'
    assert '{:E}'.format(nan) == 'NAN'
    assert '{:+E}'.format(nan) == 'NAN'

    assert '{:g}'.format(inf) == 'inf'
    assert '{:+g}'.format(inf) == '+inf'
    assert '{:G}'.format(inf) == 'INF'
    assert '{:+G}'.format(inf) == '+INF'

    assert '{:g}'.format(ninf) == '-inf'
    assert '{:+g}'.format(ninf) == '-inf'
    assert '{:G}'.format(ninf) == '-INF'
    assert '{:+G}'.format(ninf) == '-INF'

    assert '{:g}'.format(nan) == 'nan'
    assert '{:+g}'.format(nan) == 'nan'
    assert '{:G}'.format(nan) == 'NAN'
    assert '{:+G}'.format(nan) == 'NAN'
