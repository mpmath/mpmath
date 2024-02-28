from mpmath import inf, matrix, nstr


A1 = matrix([])
A2 = matrix([[]])
A3 = matrix(2)
A4 = matrix([1, 2, 3])


def test_nstr():
    m = matrix([[0.75, 0.190940654, -0.0299195971],
                [0.190940654, 0.65625, 0.205663228],
                [-0.0299195971, 0.205663228, 0.64453125e-20]])
    assert nstr(m, 4, min_fixed=-inf) == \
    '''[    0.75  0.1909                    -0.02992]
[  0.1909  0.6563                      0.2057]
[-0.02992  0.2057  0.000000000000000000006445]'''
    assert nstr(m, 4) == \
    '''[    0.75  0.1909   -0.02992]
[  0.1909  0.6563     0.2057]
[-0.02992  0.2057  6.445e-21]'''

def test_matrix_repr():
    assert repr(A1) == \
    '''matrix(
[])'''
    assert repr(A2) == \
    '''matrix(
[[]])'''
    assert repr(A3) == \
    '''matrix(
[['0.0', '0.0'],
 ['0.0', '0.0']])'''
    assert repr(A4) == \
    '''matrix(
[['1.0'],
 ['2.0'],
 ['3.0']])'''

def test_matrix_str():
    assert str(A1) == ''
    assert str(A2) == '[]'
    assert str(A3) == \
    '''[0.0  0.0]
[0.0  0.0]'''
    assert str(A4) == \
'''[1.0]
[2.0]
[3.0]'''
