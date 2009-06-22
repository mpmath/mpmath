import math
from mpmath import *

def test_bessel():
    mp.dps = 15
    assert j0(1).ae(0.765197686557966551)
    assert j0(pi).ae(-0.304242177644093864)
    assert j0(1000).ae(0.0247866861524201746)
    assert j0(-25).ae(0.0962667832759581162)
    assert j1(1).ae(0.440050585744933516)
    assert j1(pi).ae(0.284615343179752757)
    assert j1(1000).ae(0.00472831190708952392)
    assert j1(-25).ae(0.125350249580289905)
    assert besselj(5,1).ae(0.000249757730211234431)
    assert besselj(5,pi).ae(0.0521411843671184747)
    assert besselj(5,1000).ae(0.00502540694523318607)
    assert besselj(5,-25).ae(0.0660079953984229934)
    assert besselj(-3,2).ae(-0.128943249474402051)
    assert besselj(-4,2).ae(0.0339957198075684341)
    assert besselj(3,3+2j).ae(0.424718794929639595942 + 0.625665327745785804812j)
    assert besselj(0.25,4).ae(-0.374760630804249715)
    assert besselj(1+2j,3+4j).ae(0.319247428741872131 - 0.669557748880365678j)
    assert bessely(0,0) == -inf
    assert bessely(1,0) == -inf
    assert bessely(2,0) == -inf
    assert bessely(-1,0) == +inf
    assert bessely(-2,0) == -inf
    assert bessely(0,0.5).ae(-0.44451873350670655715)
    assert bessely(1,0.5).ae(-1.4714723926702430692)
    assert bessely(-1,0.5).ae(1.4714723926702430692)
    assert bessely(3.5,0.5).ae(-138.86400867242488443)
    assert bessely(0,3+4j).ae(4.6047596915010138655-8.8110771408232264208j)
    assert bessely(0,j).ae(-0.26803248203398854876+1.26606587775200833560j)
    assert besseli(0,0) == 1
    assert besseli(1,0) == 0
    assert besseli(2,0) == 0
    assert besseli(-1,0) == 0
    assert besseli(-2,0) == 0
    assert besseli(0,0.5).ae(1.0634833707413235193)
    assert besseli(1,0.5).ae(0.25789430539089631636)
    assert besseli(-1,0.5).ae(0.25789430539089631636)
    assert besseli(3.5,0.5).ae(0.00068103597085793815863)
    assert besseli(0,3+4j).ae(-3.3924877882755196097-1.3239458916287264815j)
    assert besseli(0,j).ae(besselj(0,1))
    assert besselk(0,0) == inf
    assert besselk(1,0) == inf
    assert besselk(2,0) == inf
    assert besselk(-1,0) == inf
    assert besselk(-2,0) == inf
    assert besselk(0,0.5).ae(0.92441907122766586178)
    assert besselk(1,0.5).ae(1.6564411200033008937)
    assert besselk(-1,0.5).ae(1.6564411200033008937)
    assert besselk(3.5,0.5).ae(207.48418747548460607)
    assert besselk(0,3+4j).ae(-0.007239051213570155013+0.026510418350267677215j)
    assert besselk(0,j).ae(-0.13863371520405399968-1.20196971531720649914j)

def test_hankel():
    mp.dps = 15
    assert hankel1(0,0.5).ae(0.93846980724081290423-0.44451873350670655715j)
    assert hankel1(1,0.5).ae(0.2422684576748738864-1.4714723926702430692j)
    assert hankel1(-1,0.5).ae(-0.2422684576748738864+1.4714723926702430692j)
    assert hankel1(1.5,0.5).ae(0.0917016996256513026-2.5214655504213378514j)
    assert hankel1(1.5,3+4j).ae(0.0066806866476728165382-0.0036684231610839127106j)
    assert hankel2(0,0.5).ae(0.93846980724081290423+0.44451873350670655715j)
    assert hankel2(1,0.5).ae(0.2422684576748738864+1.4714723926702430692j)
    assert hankel2(-1,0.5).ae(-0.2422684576748738864-1.4714723926702430692j)
    assert hankel2(1.5,0.5).ae(0.0917016996256513026+2.5214655504213378514j)
    assert hankel2(1.5,3+4j).ae(14.783528526098567526-7.397390270853446512j)

def test_struve():
    mp.dps = 15
    assert struveh(2,3).ae(0.74238666967748318564)
    assert struveh(-2.5,3).ae(0.41271003220971599344)
    assert struvel(2,3).ae(1.7476573277362782744)
    assert struvel(-2.5,3).ae(1.5153394466819651377)

def test_whittaker():
    mp.dps = 15
    assert whitm(2,3,4).ae(49.753745589025246591)
    assert whitw(2,3,4).ae(14.111656223052932215)

def test_kelvin():
    mp.dps = 15
    assert ber(2,3).ae(0.80836846563726819091)
    assert ber(3,4).ae(-0.28262680167242600233)
    assert ber(-3,2).ae(-0.085611448496796363669)
    assert bei(2,3).ae(-0.89102236377977331571)
    assert bei(-3,2).ae(-0.14420994155731828415)
    assert ker(2,3).ae(0.12839126695733458928)
    assert ker(-3,2).ae(-0.29802153400559142783)
    assert ker(0.5,3).ae(-0.085662378535217097524)
    assert kei(2,3).ae(0.036804426134164634000)
    assert kei(-3,2).ae(0.88682069845786731114)
    assert kei(0.5,3).ae(0.013633041571314302948)

def test_hyper_misc():
    mp.dps = 15
    assert hyp0f1(1,0) == 1
    assert hyp1f1(1,2,0) == 1
    assert hyp1f2(1,2,3,0) == 1
    assert hyp2f1(1,2,3,0) == 1
    assert hyp2f2(1,2,3,4,0) == 1
    assert hyp2f3(1,2,3,4,5,0) == 1
    # Degenerate case: 0F0
    assert hyper([],[],0) == 1
    assert hyper([],[],-2).ae(exp(-2))
    # Degenerate case: 1F0
    assert hyper([2],[],1.5) == 4
    #
    assert hyp2f1((1,3),(2,3),(5,6),mpf(27)/32).ae(1.6)
    assert hyp2f1((1,4),(1,2),(3,4),mpf(80)/81).ae(1.8)
    assert hyp2f1((2,3),(1,1),(3,2),(2+j)/3).ae(1.327531603558679093+0.439585080092769253j)
    mp.dps = 25
    v = mpc('1.2282306665029814734863026', '-0.1225033830118305184672133')
    assert hyper([(3,4),2+j,1],[1,5,j/3],mpf(1)/5+j/8).ae(v)
    mp.dps = 15

def test_elliptic_integrals():
    mp.dps = 15
    assert ellipk(0).ae(pi/2)
    assert ellipk(0.5).ae(gamma(0.25)**2/(4*sqrt(pi)))
    assert ellipk(1) == inf
    assert ellipk(1+0j) == inf
    assert ellipk(-1).ae('1.3110287771460599052')
    assert ellipk(-2).ae('1.1714200841467698589')
    assert isinstance(ellipk(-2), mpf)
    assert isinstance(ellipe(-2), mpf)
    assert ellipk(-50).ae('0.47103424540873331679')
    mp.dps = 30
    n1 = +fraction(99999,100000)
    n2 = +fraction(100001,100000)
    mp.dps = 15
    assert ellipk(n1).ae('7.1427724505817781901')
    assert ellipk(n2).ae(mpc('7.1427417367963090109', '-1.5707923998261688019'))
    assert ellipe(n1).ae('1.0000332138990829170')
    v = ellipe(n2)
    assert v.real.ae('0.999966786328145474069137')
    assert (v.imag*10**6).ae('7.853952181727432')
    assert ellipk(2).ae(mpc('1.3110287771460599052', '-1.3110287771460599052'))
    assert ellipk(50).ae(mpc('0.22326753950210985451', '-0.47434723226254522087'))
    assert ellipk(3+4j).ae(mpc('0.91119556380496500866', '0.63133428324134524388'))
    assert ellipk(3-4j).ae(mpc('0.91119556380496500866', '-0.63133428324134524388'))
    assert ellipk(-3+4j).ae(mpc('0.95357894880405122483', '0.23093044503746114444'))
    assert ellipk(-3-4j).ae(mpc('0.95357894880405122483', '-0.23093044503746114444'))
    assert isnan(ellipk(nan))
    assert isnan(ellipe(nan))
    assert ellipk(inf) == 0
    assert isinstance(ellipk(inf), mpc)
    assert ellipk(-inf) == 0
    assert ellipk(1+0j) == inf
    assert ellipe(0).ae(pi/2)
    assert ellipe(0.5).ae(pi**(mpf(3)/2)/gamma(0.25)**2 +gamma(0.25)**2/(8*sqrt(pi)))
    assert ellipe(1) == 1
    assert ellipe(1+0j) == 1
    assert ellipe(inf) == mpc(0,inf)
    assert ellipe(-inf) == inf
    assert ellipe(3+4j).ae(1.4995535209333469543-1.5778790079127582745j)
    assert ellipe(3-4j).ae(1.4995535209333469543+1.5778790079127582745j)
    assert ellipe(-3+4j).ae(2.5804237855343377803-0.8306096791000413778j)
    assert ellipe(-3-4j).ae(2.5804237855343377803+0.8306096791000413778j)
    assert ellipe(2).ae(0.59907011736779610372+0.59907011736779610372j)
    assert ellipe('1e-1000000000').ae(pi/2)
    assert ellipk('1e-1000000000').ae(pi/2)
    assert ellipe(-pi).ae(2.4535865983838923)
    mp.dps = 50
    assert ellipk(1/pi).ae('1.724756270009501831744438120951614673874904182624739673')
    assert ellipe(1/pi).ae('1.437129808135123030101542922290970050337425479058225712')
    assert ellipk(-10*pi).ae('0.5519067523886233967683646782286965823151896970015484512')
    assert ellipe(-10*pi).ae('5.926192483740483797854383268707108012328213431657645509')
    v = ellipk(pi)
    assert v.real.ae('0.973089521698042334840454592642137667227167622330325225')
    assert v.imag.ae('-1.156151296372835303836814390793087600271609993858798016')
    v = ellipe(pi)
    assert v.real.ae('0.4632848917264710404078033487934663562998345622611263332')
    assert v.imag.ae('1.0637961621753130852473300451583414489944099504180510966')
    mp.dps = 15

def test_exp_integrals():
    mp.dps = 15
    x = +e
    z = e + sqrt(3)*j
    assert ei(x).ae(8.21168165538361560)
    assert li(x).ae(1.89511781635593676)
    assert si(x).ae(1.82104026914756705)
    assert ci(x).ae(0.213958001340379779)
    assert shi(x).ae(4.11520706247846193)
    assert chi(x).ae(4.09647459290515367)
    assert fresnels(x).ae(0.437189718149787643)
    assert fresnelc(x).ae(0.401777759590243012)
    assert airyai(x).ae(0.0108502401568586681)
    assert airybi(x).ae(8.98245748585468627)
    assert ei(z).ae(3.72597969491314951 + 7.34213212314224421j)
    assert li(z).ae(2.28662658112562502 + 1.50427225297269364j)
    assert si(z).ae(2.48122029237669054 + 0.12684703275254834j)
    assert ci(z).ae(0.169255590269456633 - 0.892020751420780353j)
    assert shi(z).ae(1.85810366559344468 + 3.66435842914920263j)
    assert chi(z).ae(1.86787602931970484 + 3.67777369399304159j)
    assert fresnels(z/3).ae(0.034534397197008182 + 0.754859844188218737j)
    assert fresnelc(z/3).ae(1.261581645990027372 + 0.417949198775061893j)
    assert airyai(z).ae(-0.0162552579839056062 - 0.0018045715700210556j)
    assert airybi(z).ae(-4.98856113282883371 + 2.08558537872180623j)
    assert li(0) == 0.0
    assert li(1) == -inf
    assert li(inf) == inf
    assert isinstance(li(0.7), mpf)
    assert si(inf).ae(pi/2)
    assert si(-inf).ae(-pi/2)
    assert ci(inf) == 0
    assert ci(0) == -inf
    assert isinstance(ei(-0.7), mpf)
    assert ei(inf) == inf
    assert ei(-inf) == -0.0
    assert airyai(inf) == 0
    assert airybi(inf) == inf
    assert airyai(-inf) == 0
    assert airybi(-inf) == 0
    assert fresnels(inf) == 0.5
    assert fresnelc(inf) == 0.5
    assert fresnels(-inf) == -0.5
    assert fresnelc(-inf) == -0.5
    assert shi(0) == 0
    assert shi(inf) == inf
    assert shi(-inf) == -inf
    assert chi(0) == -inf
    assert chi(inf) == inf
    assert ei(0) == -inf
    assert ei(20+70j).ae(6.1041351911152984397e6 - 2.7324109310519928872e6j)
    assert expint(0,0) == inf
    assert expint(0,1).ae(1/e)
    assert expint(0,1.5).ae(2/exp(1.5)/3)
    assert expint(1,1).ae(-ei(-1))
    assert expint(2,0).ae(1)
    assert expint(3,0).ae(1/2.)
    assert expint(4,0).ae(1/3.)
    assert expint(-2, 0.5).ae(26/sqrt(e))
    assert expint(-1,-1) == 0
    assert expint(-2,-1).ae(-e)
    assert expint(5.5, 0).ae(2/9.)
    assert expint(2.00000001,0).ae(100000000./100000001)
    assert expint(2+3j,4-j).ae(0.0023461179581675065414+0.0020395540604713669262j)
    assert expint('1.01', '1e-1000').ae(99.9999999899412802)
    assert expint('1.000000000001', 3.5).ae(0.00697013985754701819446)
    assert expint(2,3).ae(3*ei(-3)+exp(-3))
    assert (expint(10,20)*10**10).ae(0.694439055541231353)
    assert expint(3,inf) == 0
    assert expint(3.2,inf) == 0
    assert expint(3.2+2j,inf) == 0
    assert expint(1,3j).ae(-0.11962978600800032763+0.27785620120457163717j)
    assert expint(1,3).ae(0.013048381094197037413)
    assert expint(1,-3).ae(-ei(3)-pi*j)
    assert expint(3) == expint(1,3)
    assert expint(1,-20).ae(-25615652.66405658882 - 3.1415926535897932385j)
    # tests for the asymptotic expansion
    # values checked with Mathematica ExpIntegralEi
    mp.dps = 50
    r = ei(20000)
    s = '3.8781962825045010930273870085501819470698476975019e+8681'
    assert str(r) == s
    r = ei(-200)
    s = '-6.8852261063076355977108174824557929738368086933303e-90'
    assert str(r) == s
    r =ei(20000 + 10*j)
    sre = '-3.255138234032069402493850638874410725961401274106e+8681'
    sim = '-2.1081929993474403520785942429469187647767369645423e+8681'
    assert str(r.real) == sre and str(r.imag) == sim
    mp.dps = 15
    # More asymptotic expansions
    assert chi(-10**6+100j).ae('1.3077239389562548386e+434288 + 7.6808956999707408158e+434287j')
    assert shi(-10**6+100j).ae('-1.3077239389562548386e+434288 - 7.6808956999707408158e+434287j')

def test_trig_integrals():
    mp.dps = 30
    assert si(mpf(1)/1000000).ae('0.000000999999999999944444444444446111')
    assert ci(mpf(1)/1000000).ae('-13.2382948930629912435014366276')
    assert si(10**10).ae('1.5707963267075846569685111517747537')
    assert ci(10**10).ae('-4.87506025174822653785729773959e-11')
    assert si(10**100).ae(pi/2)
    assert (ci(10**100)*10**100).ae('-0.372376123661276688262086695553')
    assert si(-3) == -si(3)
    assert ci(-3).ae(ci(3) + pi*j)

def test_airy():
    mp.dps = 15
    assert (airyai(10)*10**10).ae(1.1047532552898687)
    assert (airybi(10)/10**9).ae(0.45564115354822515)
    assert (airyai(1000)*10**9158).ae(9.306933063179556004)
    assert (airybi(1000)/10**9154).ae(5.4077118391949465477)
    assert airyai(-1000).ae(0.055971895773019918842)
    assert airybi(-1000).ae(-0.083264574117080633012)
    assert (airyai(100+100j)*10**188).ae(2.9099582462207032076 + 2.353013591706178756j)
    assert (airybi(100+100j)/10**185).ae(1.7086751714463652039 - 3.1416590020830804578j)

def test_hyper_0f1():
    mp.dps = 15
    v = 8.63911136507950465
    assert hyper([],[(1,3)],1.5).ae(v)
    assert hyper([],[1/3.],1.5).ae(v)
    assert hyp0f1(1/3.,1.5).ae(v)
    assert hyp0f1((1,3),1.5).ae(v)
    # Asymptotic expansion
    assert hyp0f1(3,1e9).ae('4.9679055380347771271e+27455')
    assert hyp0f1(3,1e9j).ae('-2.1222788784457702157e+19410 + 5.0840597555401854116e+19410j')

def test_hyper_1f1():
    mp.dps = 15
    v = 1.2917526488617656673
    assert hyper([(1,2)],[(3,2)],0.7).ae(v)
    assert hyper([(1,2)],[(3,2)],0.7+0j).ae(v)
    assert hyper([0.5],[(3,2)],0.7).ae(v)
    assert hyper([0.5],[1.5],0.7).ae(v)
    assert hyper([0.5],[(3,2)],0.7+0j).ae(v)
    assert hyper([0.5],[1.5],0.7+0j).ae(v)
    assert hyper([(1,2)],[1.5+0j],0.7).ae(v)
    assert hyper([0.5+0j],[1.5],0.7).ae(v)
    assert hyper([0.5+0j],[1.5+0j],0.7+0j).ae(v)
    assert hyp1f1(0.5,1.5,0.7).ae(v)
    assert hyp1f1((1,2),1.5,0.7).ae(v)
    # Asymptotic expansion
    assert hyp1f1(2,3,1e10).ae('2.1555012157015796988e+4342944809')
    assert (hyp1f1(2,3,1e10j)*10**10).ae(-0.97501205020039745852 - 1.7462392454512132074j)

def test_hyper_2f1():
    mp.dps = 15
    v = 1.0652207633823291032
    assert hyper([(1,2), (3,4)], [2], 0.3).ae(v)
    assert hyper([(1,2), 0.75], [2], 0.3).ae(v)
    assert hyper([0.5, 0.75], [2.0], 0.3).ae(v)
    assert hyper([0.5, 0.75], [2.0], 0.3+0j).ae(v)
    assert hyper([0.5+0j, (3,4)], [2.0], 0.3+0j).ae(v)
    assert hyper([0.5+0j, (3,4)], [2.0], 0.3).ae(v)
    assert hyper([0.5, (3,4)], [2.0+0j], 0.3).ae(v)
    assert hyper([0.5+0j, 0.75+0j], [2.0+0j], 0.3+0j).ae(v)
    v = 1.09234681096223231717 + 0.18104859169479360380j
    assert hyper([(1,2),0.75+j], [2], 0.5).ae(v)
    assert hyper([0.5,0.75+j], [2.0], 0.5).ae(v)
    assert hyper([0.5,0.75+j], [2.0], 0.5+0j).ae(v)
    assert hyper([0.5,0.75+j], [2.0+0j], 0.5+0j).ae(v)
    v = 0.9625 - 0.125j
    assert hyper([(3,2),-1],[4], 0.1+j/3).ae(v)
    assert hyper([1.5,-1.0],[4], 0.1+j/3).ae(v)
    assert hyper([1.5,-1.0],[4+0j], 0.1+j/3).ae(v)
    assert hyper([1.5+0j,-1.0+0j],[4+0j], 0.1+j/3).ae(v)
    v = 1.02111069501693445001 - 0.50402252613466859521j
    assert hyper([(2,10),(3,10)],[(4,10)],1.5).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],1.5).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],1.5+0j).ae(v)
    v = 0.76922501362865848528 + 0.32640579593235886194j
    assert hyper([(2,10),(3,10)],[(4,10)],4+2j).ae(v)
    assert hyper([0.2,(3,10)],[0.4+0j],4+2j).ae(v)
    assert hyper([0.2,(3,10)],[(4,10)],4+2j).ae(v)

def test_hyper_2f1_hard():
    mp.dps = 15
    # Singular cases
    assert hyp2f1(2,-1,-1,3) == 7
    assert hyp2f1(2,-2,-2,3) == 34
    assert hyp2f1(2,-2,-3,3) == 14
    assert hyp2f1(2,-3,-2,3) == inf
    assert hyp2f1(2,-1.5,-1.5,3) == 0.25
    assert hyp2f1(1,2,3,0) == 1
    assert hyp2f1(0,1,0,0) == 1
    assert hyp2f1(0,0,0,0) == 1
    assert isnan(hyp2f1(1,1,0,0))
    assert hyp2f1(2,-1,-5, 0.25+0.25j).ae(1.1+0.1j)
    assert hyp2f1(2,-5,-5, 0.25+0.25j).ae(163./128 + 125./128*j)
    assert hyp2f1(0.7235, -1, -5, 0.3).ae(1.04341)
    assert hyp2f1(0.7235, -5, -5, 0.3).ae(1.2939225017815903812)
    assert hyp2f1(-1,-2,4,1) == 1.5
    assert hyp2f1(1,2,-3,1) == inf
    assert hyp2f1(-2,-2,1,1) == 6
    assert hyp2f1(1,-2,-4,1).ae(5./3)
    assert hyp2f1(0,-6,-4,1) == 1
    assert hyp2f1(0,-3,-4,1) == 1
    assert hyp2f1(0,0,0,1) == 1
    assert hyp2f1(1,0,0,1) == 1
    assert hyp2f1(1,1,0,1) == inf
    assert hyp2f1(1,-6,-4,1) == inf
    assert hyp2f1(-7.2,-0.5,-4.5,1) == 0
    assert hyp2f1(-7.2,-1,-2,1).ae(-2.6)
    assert hyp2f1(1,-0.5,-4.5, 1) == inf
    assert hyp2f1(1,0.5,-4.5, 1) == -inf
    # Check evaluation on / close to unit circle
    z = exp(j*pi/3)
    w = (nthroot(2,3)+1)*exp(j*pi/12)/nthroot(3,4)**3
    assert hyp2f1('1/2','1/6','1/3', z).ae(w)
    assert hyp2f1('1/2','1/6','1/3', z.conjugate()).ae(w.conjugate())
    assert hyp2f1(0.25, (1,3), 2, '0.999').ae(1.06826449496030635)
    assert hyp2f1(0.25, (1,3), 2, '1.001').ae(1.06867299254830309446-0.00001446586793975874j)
    assert hyp2f1(0.25, (1,3), 2, -1).ae(0.96656584492524351673)
    assert hyp2f1(0.25, (1,3), 2, j).ae(0.99041766248982072266+0.03777135604180735522j)
    assert hyp2f1(2,3,5,'0.99').ae(27.699347904322690602)
    assert hyp2f1((3,2),-0.5,3,'0.99').ae(0.68403036843911661388)
    assert hyp2f1(2,3,5,1j).ae(0.37290667145974386127+0.59210004902748285917j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), 0.95*exp(j*k)) for k in range(1,15)]).ae(52.851400204289452922+6.244285013912953225j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), 1.05*exp(j*k)) for k in range(1,15)]).ae(54.506013786220655330-3.000118813413217097j)
    assert fsum([hyp2f1((7,10),(2,3),(-1,2), exp(j*k)) for k in range(1,15)]).ae(55.792077935955314887+1.731986485778500241j)
    assert hyp2f1(2,2.5,-3.25,0.999).ae(218373932801217082543180041.33)
    # Branches
    assert hyp2f1(1,1,2,1.01).ae(4.5595744415723676911-3.1104877758314784539j)
    assert hyp2f1(1,1,2,1.01+0.1j).ae(2.4149427480552782484+1.4148224796836938829j)
    assert hyp2f1(1,1,2,3+4j).ae(0.14576709331407297807+0.48379185417980360773j)
    assert hyp2f1(1,1,2,4).ae(-0.27465307216702742285 - 0.78539816339744830962j)
    assert hyp2f1(1,1,2,-4).ae(0.40235947810852509365)

def test_hyper_u():
    mp.dps = 15
    assert hyperu(2,-3,0).ae(0.05)
    assert hyperu(2,-3.5,0).ae(4./99)
    assert hyperu(2,0,0) == 0.5
    assert hyperu(-5,1,0) == -120
    assert hyperu(-5,2,0) == inf
    assert hyperu(-5,-2,0) == 0
    assert hyperu(7,7,3).ae(0.00014681269365593503986)  #exp(3)*gammainc(-6,3)
    assert hyperu(2,-3,4).ae(0.011836478100271995559)
    assert hyperu(3,4,5).ae(1./125)
    assert hyperu(2,3,0.0625) == 256
    assert hyperu(-1,2,0.25+0.5j) == -1.75+0.5j
    assert hyperu(0.5,1.5,7.25).ae(2/sqrt(29))
    assert hyperu(2,6,pi).ae(0.55804439825913399130)
    assert (hyperu((3,2),8,100+201j)*10**4).ae(-0.3797318333856738798 - 2.9974928453561707782j)
    assert (hyperu((5,2),(-1,2),-5000)*10**10).ae(-5.6681877926881664678j)
    # XXX: fails because of undetected cancellation in low level series code
    # Alternatively: could use asymptotic series here, if convergence test
    # tweaked back to recognize this one
    #assert (hyperu((5,2),(-1,2),-500)*10**7).ae(-1.82526906001593252847j)

def test_hyper_2f0():
    mp.dps = 15
    assert hyper([1,2],[],3) == hyp2f0(1,2,3)
    assert hyp2f0(2,3,7).ae(0.0116108068639728714668 - 0.0073727413865865802130j)
    assert hyp2f0(2,3,0) == 1
    assert hyp2f0(0,0,0) == 1
    assert hyp2f0(-1,-1,1).ae(2)
    assert hyp2f0(-4,1,1.5).ae(62.5)
    assert hyp2f0(-4,1,50).ae(147029801)
    assert hyp2f0(-4,1,0.0001).ae(0.99960011997600240000)
    assert hyp2f0(0.5,0.25,0.001).ae(1.0001251174078538115)
    assert hyp2f0(0.5,0.25,3+4j).ae(0.85548875824755163518 + 0.21636041283392292973j)
    # Important: cancellation check
    assert hyp2f0((1,6),(5,6),-0.02371708245126284498).ae(0.996785723120804309)
    # Should be exact; polynomial case
    assert hyp2f0(-2,1,0.5+0.5j) == 0
    assert hyp2f0(1,-2,0.5+0.5j) == 0
    # There used to be a bug in thresholds that made one of the following hang
    for d in [15, 50, 80]:
        mp.dps = d
        assert hyp2f0(1.5, 0.5, 0.009).ae('1.006867007239309717945323585695344927904000945829843527398772456281301440034218290443367270629519483 + 1.238277162240704919639384945859073461954721356062919829456053965502443570466701567100438048602352623e-46j')

def test_hyper_1f2():
    mp.dps = 15
    assert hyper([1],[2,3],4) == hyp1f2(1,2,3,4)
    a1,b1,b2 = (1,10),(2,3),1./16
    assert hyp1f2(a1,b1,b2,10).ae(298.7482725554557568)
    assert hyp1f2(a1,b1,b2,100).ae(224128961.48602947604)
    assert hyp1f2(a1,b1,b2,1000).ae(1.1669528298622675109e+27)
    assert hyp1f2(a1,b1,b2,10000).ae(2.4780514622487212192e+86)
    assert hyp1f2(a1,b1,b2,100000).ae(1.3885391458871523997e+274)
    assert hyp1f2(a1,b1,b2,1000000).ae('9.8851796978960318255e+867')
    assert hyp1f2(a1,b1,b2,10**7).ae('1.1505659189516303646e+2746')
    assert hyp1f2(a1,b1,b2,10**8).ae('1.4672005404314334081e+8685')
    assert hyp1f2(a1,b1,b2,10**20).ae('3.6888217332150976493e+8685889636')
    assert hyp1f2(a1,b1,b2,10*j).ae(-16.163252524618572878 - 44.321567896480184312j)
    assert hyp1f2(a1,b1,b2,100*j).ae(61938.155294517848171 + 637349.45215942348739j)
    assert hyp1f2(a1,b1,b2,1000*j).ae(8455057657257695958.7 + 6261969266997571510.6j)
    assert hyp1f2(a1,b1,b2,10000*j).ae(-8.9771211184008593089e+60 + 4.6550528111731631456e+59j)
    assert hyp1f2(a1,b1,b2,100000*j).ae(2.6398091437239324225e+193 + 4.1658080666870618332e+193j)
    assert hyp1f2(a1,b1,b2,1000000*j).ae('3.5999042951925965458e+613 + 1.5026014707128947992e+613j')
    assert hyp1f2(a1,b1,b2,10**7*j).ae('-8.3208715051623234801e+1939 - 3.6752883490851869429e+1941j')
    assert hyp1f2(a1,b1,b2,10**8*j).ae('2.0724195707891484454e+6140 - 1.3276619482724266387e+6141j')
    assert hyp1f2(a1,b1,b2,10**20*j).ae('-1.1734497974795488504e+6141851462 + 1.1498106965385471542e+6141851462j')

def test_hyper_2f3():
    mp.dps = 15
    assert hyper([1,2],[3,4,5],6) == hyp2f3(1,2,3,4,5,6)
    a1,a2,b1,b2,b3 = (1,10),(2,3),(3,10), 2, 1./16
    # Check asymptotic expansion
    assert hyp2f3(a1,a2,b1,b2,b3,10).ae(128.98207160698659976)
    assert hyp2f3(a1,a2,b1,b2,b3,1000).ae(6.6309632883131273141e25)
    assert hyp2f3(a1,a2,b1,b2,b3,10000).ae(4.6863639362713340539e84)
    assert hyp2f3(a1,a2,b1,b2,b3,100000).ae(8.6632451236103084119e271)
    assert hyp2f3(a1,a2,b1,b2,b3,10**6).ae('2.0291718386574980641e865')
    assert hyp2f3(a1,a2,b1,b2,b3,10**7).ae('7.7639836665710030977e2742')
    assert hyp2f3(a1,a2,b1,b2,b3,10**8).ae('3.2537462584071268759e8681')
    assert hyp2f3(a1,a2,b1,b2,b3,10**20).ae('1.2966030542911614163e+8685889627')
    assert hyp2f3(a1,a2,b1,b2,b3,10*j).ae(-18.551602185587547854 - 13.348031097874113552j)
    assert hyp2f3(a1,a2,b1,b2,b3,100*j).ae(78634.359124504488695 + 74459.535945281973996j)
    assert hyp2f3(a1,a2,b1,b2,b3,1000*j).ae(597682550276527901.59 - 65136194809352613.078j)
    assert hyp2f3(a1,a2,b1,b2,b3,10000*j).ae(-1.1779696326238582496e+59 + 1.2297607505213133872e+59j)
    assert hyp2f3(a1,a2,b1,b2,b3,100000*j).ae(2.9844228969804380301e+191 + 7.5587163231490273296e+190j)
    assert hyp2f3(a1,a2,b1,b2,b3,1000000*j).ae('7.4859161049322370311e+610 - 2.8467477015940090189e+610j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**7*j).ae('-1.7477645579418800826e+1938 - 1.7606522995808116405e+1938j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**8*j).ae('-1.6932731942958401784e+6137 - 2.4521909113114629368e+6137j')
    assert hyp2f3(a1,a2,b1,b2,b3,10**20*j).ae('-2.0988815677627225449e+6141851451 + 5.7708223542739208681e+6141851452j')

def test_hyper_2f2():
    mp.dps = 15
    assert hyper([1,2],[3,4],5) == hyp2f2(1,2,3,4,5)
    a1,a2,b1,b2 = (3,10),4,(1,2),1./16
    assert hyp2f2(a1,a2,b1,b2,10).ae(448225936.3377556696)
    assert hyp2f2(a1,a2,b1,b2,10000).ae('1.2012553712966636711e+4358')
    assert hyp2f2(a1,a2,b1,b2,-20000).ae(-0.04182343755661214626)
    assert hyp2f2(a1,a2,b1,b2,10**20).ae('1.1148680024303263661e+43429448190325182840')

def test_orthpoly():
    mp.dps = 15
    assert jacobi(-4,2,3,0.7).ae(22800./4913)
    assert jacobi(3,2,4,5.5) == 4133.125
    assert jacobi(1.5,5/6.,4,0).ae(-1.0851951434075508417)
    assert legendre(5, 7) == 129367
    assert legendre(0.5,0).ae(0.53935260118837935667)
    assert legendre(-1,-1) == 1
    assert legendre(0,-1) == 1
    assert legendre(0, 1) == 1
    assert legendre(1, -1) == -1
    assert legendre(7, 1) == 1
    assert legendre(7, -1) == -1
    assert legendre(8,1.5).ae(15457523./32768)
    assert legendre(j,-j).ae(2.4448182735671431011 + 0.6928881737669934843j)
    assert chebyu(5,1) == 6
    assert chebyt(3,2) == 26
    assert legendre(3.5,-1) == inf
    assert legendre(4.5,-1) == -inf
    assert legendre(3.5+1j,-1) == mpc(inf,inf)
    assert legendre(4.5+1j,-1) == mpc(-inf,-inf)

def test_agm():
    mp.dps = 15
    assert agm(0,0) == 0
    assert agm(0,1) == 0
    assert agm(1,1) == 1
    assert agm(7,7) == 7
    assert agm(j,j) == j
    assert (1/agm(1,sqrt(2))).ae(0.834626841674073186)
    assert agm(1,2).ae(1.4567910310469068692)
    assert agm(1,3).ae(1.8636167832448965424)
    assert agm(1,j).ae(0.599070117367796104+0.599070117367796104j)
    assert agm(2) == agm(1,2)
    assert agm(-3,4).ae(0.63468509766550907+1.3443087080896272j)

def test_incomplete_gamma():
    mp.dps = 15
    assert gammainc(2,5).ae(6*exp(-5))
    assert gammainc(2,0,5).ae(1-6*exp(-5))
    assert gammainc(2,3,5).ae(-6*exp(-5)+4*exp(-3))
    assert gammainc(-2.5,-0.5).ae(-0.9453087204829418812-5.3164237738936178621j)
    assert gammainc(0,2,4).ae(0.045121158298212213088)
    assert gammainc(0,3).ae(0.013048381094197037413)
    assert gammainc(0,2+j,1-j).ae(0.00910653685850304839-0.22378752918074432574j)
    assert gammainc(0,1-j).ae(0.00028162445198141833+0.17932453503935894015j)
    assert gammainc(3,4,5,True).ae(0.11345128607046320253)
    assert gammainc(3.5,0,inf).ae(gamma(3.5))

def test_incomplete_beta():
    mp.dps = 15
    assert betainc(-2,-3,0.5,0.75).ae(63.4305673311255413583969)
    assert betainc(4.5,0.5+2j,2.5,6).ae(0.2628801146130621387903065 + 0.5162565234467020592855378j)
    assert betainc(4,5,0,6).ae(90747.77142857142857142857)

def test_erf():
    mp.dps = 15
    assert erf(0) == 0
    assert erf(1).ae(0.84270079294971486934)
    assert erf(3+4j).ae(-120.186991395079444098 - 27.750337293623902498j)
    assert erf(-4-3j).ae(-0.99991066178539168236 + 0.00004972026054496604j)
    assert erf(pi).ae(0.99999112385363235839)
    assert erf(1j).ae(1.6504257587975428760j)
    assert erf(-1j).ae(-1.6504257587975428760j)
    assert isinstance(erf(1), mpf)
    assert isinstance(erf(-1), mpf)
    assert isinstance(erf(0), mpf)
    assert isinstance(erf(0j), mpc)
    assert erf(inf) == 1
    assert erf(-inf) == -1
    assert erfi(0) == 0
    assert erfi(1/pi).ae(0.371682698493894314)
    assert erfi(inf) == inf
    assert erfi(-inf) == -inf
    assert erf(1+0j) == erf(1)
    assert erfc(1+0j) == erfc(1)
    assert erf(0.2+0.5j).ae(1 - erfc(0.2+0.5j))
    assert erfc(0) == 1
    assert erfc(1).ae(1-erf(1))
    assert erfc(-1).ae(1-erf(-1))
    assert erfc(1/pi).ae(1-erf(1/pi))
    assert erfc(-10) == 2
    assert erfc(-1000000) == 2
    assert erfc(-inf) == 2
    assert erfc(inf) == 0
    assert isnan(erfc(nan))
    assert (erfc(10**4)*mpf(10)**43429453).ae('3.63998738656420')
    mp.dps = 50
    # This one does not use the asymptotic series
    assert (erfc(10)*10**45).ae('2.0884875837625447570007862949577886115608181193212')
    # This one does
    assert (erfc(50)*10**1088).ae('2.0709207788416560484484478751657887929322509209954')
    mp.dps = 15
    assert str(erfc(10**50)) == '3.66744826532555e-4342944819032518276511289189166050822943970058036665661144537831658646492088707747292249493384317534'
    assert erfinv(0) == 0
    assert erfinv(0.5).ae(0.47693627620446987338)
    assert erfinv(-0.5).ae(-0.47693627620446987338)
    assert erfinv(1) == inf
    assert erfinv(-1) == -inf
    assert erf(erfinv(0.95)).ae(0.95)
    assert erf(erfinv(0.999999999995)).ae(0.999999999995)
    assert erf(erfinv(-0.999999999995)).ae(-0.999999999995)
    mp.dps = 50
    assert erf(erfinv('0.99999999999999999999999999999995')).ae('0.99999999999999999999999999999995')
    assert erf(erfinv('0.999999999999999999999999999999995')).ae('0.999999999999999999999999999999995')
    assert erf(erfinv('-0.999999999999999999999999999999995')).ae('-0.999999999999999999999999999999995')
    mp.dps = 15
    # Complex asymptotic expansions
    v = erfc(50j)
    assert v.real == 1
    assert v.imag.ae('-6.1481820666053078736e+1083')
    assert erfc(-100+5j).ae(2)
    assert (erfc(100+5j)*10**4335).ae(2.3973567853824133572 - 3.9339259530609420597j)
    assert erfc(100+100j).ae(0.00065234366376857698698 - 0.0039357263629214118437j)

def test_pdf():
    mp.dps = 15
    assert npdf(-inf) == 0
    assert npdf(inf) == 0
    assert npdf(5,0,2).ae(5+4,4,2)
    assert quadts(lambda x: npdf(x,-0.5,0.8), [-inf, inf]) == 1
    assert ncdf(0) == 0.5
    assert ncdf(3,3) == 0.5
    assert ncdf(-inf) == 0
    assert ncdf(inf) == 1
    assert ncdf(10) == 1
    # Verify that this is computed accurately
    assert (ncdf(-10)*10**24).ae(7.619853024160526)

def test_lambertw():
    mp.dps = 15
    assert lambertw(0) == 0
    assert lambertw(0+0j) == 0
    assert lambertw(inf) == inf
    assert isnan(lambertw(nan))
    assert lambertw(inf,1).real == inf
    assert lambertw(inf,1).imag.ae(2*pi)
    assert lambertw(-inf,1).real == inf
    assert lambertw(-inf,1).imag.ae(3*pi)
    assert lambertw(0,-1) == -inf
    assert lambertw(0,1) == -inf
    assert lambertw(0,3) == -inf
    assert lambertw(e).ae(1)
    assert lambertw(1).ae(0.567143290409783873)
    assert lambertw(-pi/2).ae(j*pi/2)
    assert lambertw(-log(2)/2).ae(-log(2))
    assert lambertw(0.25).ae(0.203888354702240164)
    assert lambertw(-0.25).ae(-0.357402956181388903)
    assert lambertw(-1./10000,0).ae(-0.000100010001500266719)
    assert lambertw(-0.25,-1).ae(-2.15329236411034965)
    assert lambertw(0.25,-1).ae(-3.00899800997004620-4.07652978899159763j)
    assert lambertw(-0.25,-1).ae(-2.15329236411034965)
    assert lambertw(0.25,1).ae(-3.00899800997004620+4.07652978899159763j)
    assert lambertw(-0.25,1).ae(-3.48973228422959210+7.41405453009603664j)
    assert lambertw(-4).ae(0.67881197132094523+1.91195078174339937j)
    assert lambertw(-4,1).ae(-0.66743107129800988+7.76827456802783084j)
    assert lambertw(-4,-1).ae(0.67881197132094523-1.91195078174339937j)
    assert lambertw(1000).ae(5.24960285240159623)
    assert lambertw(1000,1).ae(4.91492239981054535+5.44652615979447070j)
    assert lambertw(1000,-1).ae(4.91492239981054535-5.44652615979447070j)
    assert lambertw(1000,5).ae(3.5010625305312892+29.9614548941181328j)
    assert lambertw(3+4j).ae(1.281561806123775878+0.533095222020971071j)
    assert lambertw(-0.4+0.4j).ae(-0.10396515323290657+0.61899273315171632j)
    assert lambertw(3+4j,1).ae(-0.11691092896595324+5.61888039871282334j)
    assert lambertw(3+4j,-1).ae(0.25856740686699742-3.85211668616143559j)
    assert lambertw(-0.5,-1).ae(-0.794023632344689368-0.770111750510379110j)
    assert lambertw(-1./10000,1).ae(-11.82350837248724344+6.80546081842002101j)
    assert lambertw(-1./10000,-1).ae(-11.6671145325663544)
    assert lambertw(-1./10000,-2).ae(-11.82350837248724344-6.80546081842002101j)
    assert lambertw(-1./100000,4).ae(-14.9186890769540539+26.1856750178782046j)
    assert lambertw(-1./100000,5).ae(-15.0931437726379218666+32.5525721210262290086j)
    assert lambertw((2+j)/10).ae(0.173704503762911669+0.071781336752835511j)
    assert lambertw((2+j)/10,1).ae(-3.21746028349820063+4.56175438896292539j)
    assert lambertw((2+j)/10,-1).ae(-3.03781405002993088-3.53946629633505737j)
    assert lambertw((2+j)/10,4).ae(-4.6878509692773249+23.8313630697683291j)
    assert lambertw(-(2+j)/10).ae(-0.226933772515757933-0.164986470020154580j)
    assert lambertw(-(2+j)/10,1).ae(-2.43569517046110001+0.76974067544756289j)
    assert lambertw(-(2+j)/10,-1).ae(-3.54858738151989450-6.91627921869943589j)
    assert lambertw(-(2+j)/10,4).ae(-4.5500846928118151+20.6672982215434637j)
    mp.dps = 50
    assert lambertw(pi).ae('1.073658194796149172092178407024821347547745350410314531')
    mp.dps = 15
