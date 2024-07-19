import mpmath as mp

try:
    import _testcapi
except ImportError:
    _testcapi = None

INF = float("inf")
NAN = float("nan")


# locate files with float format test values
format_testfile = 'format_mpf_testcases.txt'
format_testfile_2 = 'format_mpf_testcases_2.txt'


def test_mpf_fmt_cpython():
    '''
    These tests assure that mpf.__format__ yields the same result as regular
    float.__format__, when dps is default.
    '''
    with open(format_testfile, encoding="utf-8") as testfile:
        for line in testfile:
            if line.startswith('--') or line.startswith('%r'):
                continue
            line = line.strip()
            if not line:
                continue

            lhs, rhs = map(str.strip, line.split('->'))
            fmt, arg = lhs.split()
            f = mp.mpf(arg)

            if fmt != '%r':
                fmt2 = fmt[1:]
                assert format(f, fmt2) == rhs
                # Negative 0 is not covered by mpmath

                if float(rhs) != 0:
                    assert format(-f, fmt2) == '-' + rhs


def test_mpf_fmt():
    '''
    These tests are specific to mpf.
    '''

    with mp.workdps(1000):
        with open(format_testfile_2, encoding="utf-8") as testfile:
            for line in testfile:
                if line.startswith('--'):
                    continue
                line = line.strip()
                if not line:
                    continue

                lhs, rhs = map(str.strip, line.split('->'))
                fmt, arg = lhs.split()

                f = mp.mpf(arg)
                assert format(f, fmt) == rhs


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
