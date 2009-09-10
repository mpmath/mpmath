from usertools import monitor, timing

from ctx_fp import FPContext
from ctx_mp import MPContext

fp = FPContext()
mp = MPContext()

bernfrac = mp.bernfrac

jdn = mp.jdn
jsn = mp.jsn
jcn = mp.jcn
jtheta = mp.jtheta

plot = mp.plot
cplot = mp.cplot

odefun = mp.odefun

jacobian = mp.jacobian
findroot = mp.findroot
multiplicity = mp.multiplicity

isinf = mp.isinf
isnan = mp.isnan
nan = mp.nan
rand = mp.rand

fraction = mp.fraction

linspace = mp.linspace
arange = mp.arange

mpmathify = convert = mp.convert
mpc = mp.mpc

nstr = mp.nstr
nprint = mp.nprint
chop = mp.chop

fmul = mp.fmul
fprod = mp.fprod

quad = mp.quad
quadgl = mp.quadgl
quadts = mp.quadts
quadosc = mp.quadosc

pslq = mp.pslq
identify = mp.identify
findpoly = mp.findpoly

richardson = mp.richardson
shanks = mp.shanks
nsum = mp.nsum
nprod = mp.nprod
diff = mp.diff
diffs = mp.diffs
diffun = mp.diffun
differint = mp.differint
taylor = mp.taylor
pade = mp.pade
polyval = mp.polyval
polyroots = mp.polyroots
fourier = mp.fourier
fourierval = mp.fourierval
sumem = mp.sumem
chebyfit = mp.chebyfit
limit = mp.limit

matrix = mp.matrix
eye = mp.eye
diag = mp.diag
zeros = mp.zeros
ones = mp.ones
hilbert = mp.hilbert
randmatrix = mp.randmatrix
swap_row = mp.swap_row
extend = mp.extend
norm = mp.norm
mnorm = mp.mnorm

lu_solve = mp.lu_solve
lu = mp.lu
unitvector = mp.unitvector
inverse = mp.inverse
residual = mp.residual
qr_solve = mp.qr_solve
cholesky = mp.cholesky
cholesky_solve = mp.cholesky_solve
det = mp.det
cond = mp.cond
exp_pade = mp.exp_pade

mpf = mp.mpf
j = mp.j
exp = mp.exp
ln = mp.ln
im = mp.im
re = mp.re
inf = mp.inf
ninf = mp.ninf
lambertw = mp.lambertw
sign = mp.sign

eps = mp.eps
pi = mp.pi
ln2 = mp.ln2
ln10 = mp.ln10
phi = mp.phi
e = mp.e
euler = mp.euler
catalan = mp.catalan
khinchin = mp.khinchin
glaisher = mp.glaisher
apery = mp.apery
degree = mp.degree
twinprime = mp.twinprime
mertens = mp.mertens

ldexp = mp.ldexp
frexp = mp.frexp

fsum = mp.fsum
fdot = mp.fdot

sqrt = mp.sqrt
cbrt = mp.cbrt
exp = mp.exp
ln = mp.ln
log = mp.log
log10 = mp.log10
power = mp.power
cos = mp.cos
sin = mp.sin
tan = mp.tan
cosh = mp.cosh
sinh = mp.sinh
tanh = mp.tanh
acos = mp.acos
asin = mp.asin
atan = mp.atan
asinh = mp.asinh
acosh = mp.acosh
atanh = mp.atanh
sec = mp.sec
csc = mp.csc
cot = mp.cot
sech = mp.sech
csch = mp.csch
coth = mp.coth
asec = mp.asec
acsc = mp.acsc
acot = mp.acot
asech = mp.asech
acsch = mp.acsch
acoth = mp.acoth
cospi = mp.cospi
sinpi = mp.sinpi
sinc = mp.sinc
sincpi = mp.sincpi
fabs = mp.fabs
re = mp.re
im = mp.im
conj = mp.conj
floor = mp.floor
ceil = mp.ceil
root = mp.root
nthroot = mp.nthroot
hypot = mp.hypot
modf = mp.modf
ldexp = mp.ldexp
frexp = mp.frexp
sign = mp.sign
arg = mp.arg
degrees = mp.degrees
radians = mp.radians
atan2 = mp.atan2
fib = mp.fib
fibonacci = mp.fibonacci
zeta = mp.zeta
altzeta = mp.altzeta
gamma = mp.gamma
factorial = mp.factorial
fac = mp.fac
fac2 = mp.fac2
beta = mp.beta
betainc = mp.betainc
psi = mp.psi
psi0 = mp.psi0
psi1 = mp.psi1
psi2 = mp.psi2
psi3 = mp.psi3
polygamma = mp.polygamma
digamma = mp.digamma
trigamma = mp.trigamma
tetragamma = mp.tetragamma
pentagamma = mp.pentagamma
harmonic = mp.harmonic
bernoulli = mp.bernoulli
bernfrac = mp.bernfrac
stieltjes = mp.stieltjes
hurwitz = mp.hurwitz
dirichlet = mp.dirichlet
bernpoly = mp.bernpoly
polylog = mp.polylog
clsin = mp.clsin
clcos = mp.clcos
gammainc = mp.gammainc
gammaprod = mp.gammaprod
binomial = mp.binomial
rf = mp.rf
ff = mp.ff
hyper = mp.hyper
hyp0f1 = mp.hyp0f1
hyp1f1 = mp.hyp1f1
hyp1f2 = mp.hyp1f2
hyp2f1 = mp.hyp2f1
hyp2f2 = mp.hyp2f2
hyp2f0 = mp.hyp2f0
hyp2f3 = mp.hyp2f3
hyperu = mp.hyperu
hypercomb = mp.hypercomb
meijerg = mp.meijerg
appellf1 = mp.appellf1
erf = mp.erf
erfc = mp.erfc
erfi = mp.erfi
erfinv = mp.erfinv
npdf = mp.npdf
ncdf = mp.ncdf
expint = mp.expint
e1 = mp.e1
ei = mp.ei
li = mp.li
ci = mp.ci
si = mp.si
chi = mp.chi
shi = mp.shi
fresnels = mp.fresnels
fresnelc = mp.fresnelc
airyai = mp.airyai
airybi = mp.airybi
ellipe = mp.ellipe
ellipk = mp.ellipk
agm = mp.agm
jacobi = mp.jacobi
chebyt = mp.chebyt
chebyu = mp.chebyu
legendre = mp.legendre
legenp = mp.legenp
legenq = mp.legenq
hermite = mp.hermite
gegenbauer = mp.gegenbauer
laguerre = mp.laguerre
besselj = mp.besselj
j0 = mp.j0
j1 = mp.j1
besseli = mp.besseli
bessely = mp.bessely
besselk = mp.besselk
hankel1 = mp.hankel1
hankel2 = mp.hankel2
struveh = mp.struveh
struvel = mp.struvel
whitm = mp.whitm
whitw = mp.whitw
ber = mp.ber
bei = mp.bei
ker = mp.ker
kei = mp.kei
coulombc = mp.coulombc
coulombf = mp.coulombf
coulombg = mp.coulombg
lambertw = mp.lambertw
barnesg = mp.barnesg
superfac = mp.superfac
hyperfac = mp.hyperfac
loggamma = mp.loggamma
siegeltheta = mp.siegeltheta
siegelz = mp.siegelz
grampoint = mp.grampoint
zetazero = mp.zetazero
riemannr = mp.riemannr
primepi = mp.primepi
primepi2 = mp.primepi2
primezeta = mp.primezeta
bell = mp.bell
polyexp = mp.polyexp
expm1 = mp.expm1
powm1 = mp.powm1
unitroots = mp.unitroots
cyclotomic = mp.cyclotomic





"""

from mp4 import *
mp.lambertw(3)
fp.lambertw(3)


import os
os.environ['MPMATH_NOGMPY'] = '1'
from mp4 import fsum, timing
fsum([2,3,4])
1/timing(fsum, [2,3,4])


"""


if __name__ == '__main__':
    #import doctest
    #doctest.testmod()
    try:
        import psyco; psyco.full()
    except ImportError:
        pass
    import sys
    from timeit import default_timer as clock
    filter = []
    for i, arg in enumerate(sys.argv):
        if '__init__.py' in arg:
            filter = sys.argv[i+1:]
            break
    import doctest
    globs = globals().copy()
    for obj in globs: #sorted(globs.keys()):
        if filter:
            if not sum([pat in obj for pat in filter]):
                continue
        print obj,
        t1 = clock()
        doctest.run_docstring_examples(globs[obj], {})
        t2 = clock()
        print round(t2-t1, 3)

