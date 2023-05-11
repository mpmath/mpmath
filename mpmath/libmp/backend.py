import os
import sys

#----------------------------------------------------------------------------#
# Support GMPY for high-speed large integer arithmetic.                      #
#                                                                            #
# To allow an external module to handle arithmetic, we need to make sure     #
# that all high-precision variables are declared of the correct type. MPZ    #
# is the constructor for the high-precision type. It defaults to Python's    #
# long type but can be assinged another type, typically gmpy.mpz.            #
#                                                                            #
# MPZ must be used for the mantissa component of an mpf and must be used     #
# for internal fixed-point operations.                                       #
#                                                                            #
# Side-effects                                                               #
# 1) "is" cannot be used to test for special values. Must use "==".          #
# 2) There are bugs in GMPY prior to v1.02 so we must use v1.03 or later.    #
#----------------------------------------------------------------------------#

# So we can import it from this module
gmpy = None
sage = None
sage_utils = None
BACKEND = 'python'
MPZ = int
HASH_MODULUS = sys.hash_info.modulus
HASH_BITS = 31 if sys.hash_info.width == 32 else 61


if 'MPMATH_NOGMPY' not in os.environ:
    try:
        import gmpy2 as gmpy
        BACKEND = 'gmpy'
        MPZ = gmpy.mpz
    except ImportError:
        pass

if ('MPMATH_NOSAGE' not in os.environ and 'SAGE_ROOT' in os.environ or
        'MPMATH_SAGE' in os.environ):
    try:
        import sage.all
        import sage.libs.mpmath.utils as _sage_utils
        sage = sage.all
        sage_utils = _sage_utils
        BACKEND = 'sage'
        MPZ = sage.Integer
    except:
        pass

if 'MPMATH_STRICT' in os.environ:
    STRICT = True
else:
    STRICT = False

MPZ_ZERO = MPZ(0)
MPZ_ONE = MPZ(1)
MPZ_TWO = MPZ(2)
MPZ_THREE = MPZ(3)
MPZ_FIVE = MPZ(5)

int_types = (int,) if BACKEND == 'python' else (int, MPZ)
