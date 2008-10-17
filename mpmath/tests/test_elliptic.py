"""
    Limited tests of the elliptic functions module.  A full suite of
    extensive testing can be found in elliptic_torture_tests.py

    Author of the first version: M.T. Taschuk

    References:

    [1] Abramowitz & Stegun. 'Handbook of Mathematical Functions, 9th Ed.', 
        (Dover duplicate of 1972 edition)
    [2] Whittaker 'A Course of Modern Analysis, 4th Ed.', 1946, 
        Cambridge Univeristy Press

"""
__version__ = '$Id:$'

import unittest
#import mpmath.mptypes      # switch to this once integrated with mpmath
import mpmath
import random
from mpmath.mptypes import (mpc, eps, j)

def mpc_ae(a, b, eps=eps):
    res = True
    res = res and a.real.ae(b.real, eps)
    res = res and a.imag.ae(b.imag, eps)
    return res


from mpmath.elliptic import *

class precisemathTests(unittest.TestCase):

    def runTest(self):
        pass

    def testCalculateNome(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        q = calculate_nome(mpmath.mpf('0'))
        self.assertEquals(mpmath.mpf('0'), mpmath.mpf('0'))

        mathematica = [ (0.1,   0.00658465), 
                        (0.3,   0.0222774),
                        (0.5,   0.0432139),
                        (0.7,   0.0746899),
                        (0.9,   0.140173),
                        (0.99,  0.262196)]
        
        for i in mathematica:
            m = mpmath.mpf(i[0])
            value = calculate_nome(m.sqrt())
            self.assertEquals(round(i[1], 6), round(value, 6))
    
    def testJacobiTheta1(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        #self.assertRaises(TypeError, jacobi_theta_1, 0, 0)
        
        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_1(z, q)
        self.assertEquals(value, mpmath.mpf('0'))

        z = mpmath.mpf('0')
        m = mpmath.pi 
        self.assertRaises(ValueError, jacobi_theta_1, z, m)
        m = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, m)

        # Mathematica value for v1(u = 0.1, q = 0.1) = 0.108958
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_1(z, m) 
        self.assertEquals(round(result, 6), 0.108958)
        self.assertTrue(isinstance(result, mpmath.mpf))
        
    def testJacobiTheta2(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')

        #self.assertRaises(TypeError, jacobi_theta_2, 0, 0)

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_2(z, q)
        self.assertEquals(value, mpmath.mpf('0'))

        # Mathematica value for v2(z = 0.1, q = 0.1) = 1.12981
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_2(z, m)     # verbosity on
        self.assertEquals(round(result, 5), 1.12981)
       
        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_2(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def testJacobiTheta3(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        one = mpmath.mpf('1')
        
        #self.assertRaises(TypeError, jacobi_theta_2, 0, 0)

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_3(z, q)
        self.assertEquals(mpmath.mpf('1'), value)
        self.assertTrue(isinstance(value, mpmath.mpf))

        # Mathematica value for v3(z = 0.1, q = 0.1) = 1.1962
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_3(z, m)
        self.assertEquals(round(result, 4), 1.1962)

        mpmath.mpf.dps = 2
        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_2(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def testJacobiTheta4(self):
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 2))
        #print >> sys.stderr, testlimit

        #self.assertRaises(TypeError, jacobi_theta_4, 0, 0)

        z = mpmath.mpf('0')
        q = mpmath.mpf('0')
        value = jacobi_theta_4(z, q)
        self.assertEquals(value, mpmath.mpf('1.0'))

        # Mathematica value for v4(z = 0.1, q = 0.1) = 0.804171
        # q = 0.1, therefore m = 0.802403, according to Mathematica
        z = mpmath.mpf('0.1')
        m = mpmath.mpf('0.802403')

        result = jacobi_theta_4(z, m)
        self.assertEquals(round(result, 6), 0.804171)
       
        z = mpmath.mpf('0')
        q = mpmath.pi / mpmath.mpf('2')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)
        q = mpmath.mpf('1')
        self.assertRaises(ValueError, jacobi_theta_1, z, q)

        z = mpmath.mpf('0.1')
        q = mpmath.mpf('0.1')
        value = jacobi_theta_4(z, q)
        self.assertTrue(isinstance(value, mpmath.mpf))

    def test_JacobiTheta_complex(self):
        mp.dps = 30
        z = mpf(1)/4 + j/8
        q = mpf(1)/3 + j/7
        #N[EllipticTheta[1, 1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('0.31618034835986160705729105731678285') + \
              mpf('0.07542013825835103435142515194358975') * j
        r = jacobi_theta(1, z, q)
        self.assertTrue(mpc_ae(r, res))
    
        #N[EllipticTheta[2, 1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('1.6530986428239765928634711417951828') + \
              mpf('0.2015344864707197230526742145361455') * j
        r = jacobi_theta(2, z, q)
        self.assertTrue(mpc_ae(r, res))
    
        #N[EllipticTheta[3, 1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('1.6520564411784228184326012700348340') + \
              mpf('0.1998129119671271328684690067401823') * j
        r = jacobi_theta(3, z, q)
        self.assertTrue(mpc_ae(r, res))
    
        #N[EllipticTheta[4, 1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('0.37619082382228348252047624089973824') - \
              mpf('0.15623022130983652972686227200681074') * j
        r = jacobi_theta(4, z, q)
        self.assertTrue(mpc_ae(r, res))
    
        # check that the two representations of the theta function
        # are consistent with each other
        m = mpf(1)/3 + j/5
        q = calculate_nome(sqrt(m))
        r1 = jacobi_theta_1(z, m)
        r2 = jacobi_theta(1, z, q)
        self.assertTrue(mpc_ae(r1, r2))
        r1 = jacobi_theta_2(z, m)
        r2 = jacobi_theta(2, z, q)
        self.assertTrue(mpc_ae(r1, r2))
        r1 = jacobi_theta_3(z, m)
        r2 = jacobi_theta(3, z, q)
        self.assertTrue(mpc_ae(r1, r2))
        r1 = jacobi_theta_4(z, m)
        r2 = jacobi_theta(4, z, q)
        self.assertTrue(mpc_ae(r1, r2))

        # check some theta function identities
        theta = jacobi_theta
        mp.dos = 100
        z = mpf(1)/4 + j/8
        q = mpf(1)/3 + j/7
        mp.dps += 10
        a = [0,0,theta(2, 0, q), theta(3, 0, q), theta(4, 0, q)]
        t = [0, theta(1, z, q), theta(2, z, q), theta(3, z, q), theta(4, z, q)]
        r = [(t[2]*a[4])**2 - (t[4]*a[2])**2 + (t[1] *a[3])**2,
            (t[3]*a[4])**2 - (t[4]*a[3])**2 + (t[1] *a[2])**2,
            (t[1]*a[4])**2 - (t[3]*a[2])**2 + (t[2] *a[3])**2,
            (t[4]*a[4])**2 - (t[3]*a[3])**2 + (t[2] *a[2])**2,
            a[2]**4 + a[4]**4 - a[3]**4]
        mp.dps -= 10
        for x in r:
            self.assertTrue(mpc_ae(x, mpc(0)))


    def testJacobiEllipticSn(self):
        """
        Test some special cases of the sn(z, q) function.

        This is an intensive test, so precision turned down during
        development.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')
        
        # trival case
        result = jacobi_elliptic_sn(zero, zero)
        self.assertEquals(result, zero)

        # Abramowitz Table 16.5
        #
        # sn(K, m) = 1; K is K(k), first complete elliptic integral
        mstring = str(random.random())
        m = mpmath.mpf(mstring)
        k = m.sqrt()

        K = mpmath.ellipk(k**2)

        equality = abs(one - jacobi_elliptic_sn(K, m))

        if equality < testlimit:
            self.assertEquals(True, True)
        else:
            print >> sys.stderr, 'Sn, sn(K, m) - 1 == 0: %e' % equality
            equality = jacobi_elliptic_sn(K, m, True)   # verbose
            self.assertEquals(False, True)

        # Abramowitz Table 16.6.1
        #
        # sn(z, 0) = sin(z), m == 0
        #
        # sn(z, 1) = tanh(z), m == 1
        #
        # It would be nice to test these, but I find that they run
        # in to numerical trouble.  I'm currently treating as a boundary
        # case for sn function.
        
        # Mathematica value for sn(z = 0.1, m = 0.1) = 0.0998169
        arg = mpmath.mpf('0.1')
        result = jacobi_elliptic_sn(arg, arg)
        self.assertEquals(round(result, 7), 0.0998169)

    def testJacobiEllipticCn(self):
        """
        Test some special cases of the cn(z, q) function.

        This is an intensive test, so precision turned down during
        development.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # Abramowitz Table 16.5
        #
        # cn(0, q) = 1

        qstring = str(random.random())
        q = mpmath.mpf(qstring)

        cn = jacobi_elliptic_cn(zero, q)
        equality = one - cn

        if equality < testlimit:
            self.assertEquals(True, True)
        else:
            print >> sys.stderr, 'Cn (~ 1): %e' % cn
            print >> sys.stderr, 'Equality (~ zero): %e' % equality
            self.assertEquals(False, True)

        # Abramowitz Table 16.5
        #
        # cn(K, q) = 0; K is K(k), first complete elliptic integral

        mstring = str(random.random())
        m = mpmath.mpf(mstring)
        k = m.sqrt()

        K = mpmath.ellipk(k**2)

        equality = jacobi_elliptic_cn(K, m)

        if equality < testlimit:
            self.assertEquals(True, True)
        else:
            print >> sys.stderr, '\n**** Cn failure ****'
            print >> sys.stderr, '\nK: %e' % K,
            print >> sys.stderr, '\tm: %f' % m, 
            print >> sys.stderr, '\tcn: %e' % equality
            equality = jacobi_elliptic_cn(K, k, True)
            self.assertEquals(False, True)

        # Abramowitz Table 16.6.2
        #
        # cn(u, 0) = cos(u), m == 0
        #
        # cn(u, 1) = sech(z), m == 1
        #
        # It would be nice to test these, but I find that they run
        # in to numerical trouble.  I'm currently treating as a boundary
        # case for cn function.

    def testJacobiEllipticDn(self):
        """
        Test some special cases of the dn(z, q) function.
        """
        mpmath.mpf.dps = 100
        mpmath.mp.dps = 100
        #mpmath.mpf.dps = 20             # testing version
        #mpmath.mp.dps = 20
        testlimit = mpmath.mpf('10')**(-1*(mpmath.mpf.dps - 4))
        #print >> sys.stderr, testlimit

        zero = mpmath.mpf('0')
        one = mpmath.mpf('1')

        # Abramowitz Table 16.5
        #
        # dn(0, q) = 1

        mstring = str(random.random())
        m = mpmath.mpf(mstring)

        dn = jacobi_elliptic_dn(zero, m)
        equality = one - dn

        if equality < testlimit:
            self.assertEquals(True, True)
        else:
            print >> sys.stderr, '\n**** Dn failure ****'
            print >> sys.stderr, '\tm: %f' % m, 
            print >> sys.stderr, '\tdn: %e' % dn,
            print >> sys.stderr, '\tEquality: %e' % equality
            equality = jacobi_elliptic_dn(zero, m, True)
            self.assertEquals(False, True)

        # Abramowitz Table 16.6.3
        #
        # dn(z, 0) = 1, m == 0

        zstring = str(random.random())
        z = mpmath.mpf(zstring)

        value = jacobi_elliptic_dn(z, zero)

        equality = value - one

        if equality < testlimit:
            self.assertEquals(True, True)
        else:
            print >> sys.stderr, 'Equality (~ zero): %e' % equality
            self.assertEquals(False, True)

    def test_JacobiElliptic_complex(self):
        mp.dps = 30
        # N[JacobiSN[1/4 + I/8, 1/3 + I/7], 35] in Mathematica
        res = mpf('0.2495674401066275492326652143537') + \
              mpf('0.12017344422863833381301051702823') * j
        u = mpf(1)/4 + j/8
        m = mpf(1)/3 + j/7
        r = jacobi_elliptic_sn(u, m)
        self.assertTrue(mpc_ae(r, res))
    
        #N[JacobiCN[1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('0.9762691700944007312693721148331') - \
              mpf('0.0307203994181623243583169154824')*j
        r = jacobi_elliptic_cn(u, m)
        #assert r.real.ae(res.real)
        #assert r.imag.ae(res.imag)
        self.assertTrue(mpc_ae(r, res))
    
        #N[JacobiDN[1/4 + I/8, 1/3 + I/7], 35]
        res = mpf('0.99639490163039577560547478589753039') - \
              mpf('0.01346296520008176393432491077244994')*j
        r = jacobi_elliptic_dn(u, m)
        self.assertTrue(mpc_ae(r, res))

def test_elliptic_functions():
    t = precisemathTests()
    t.testCalculateNome()
    t.testJacobiTheta1()
    t.testJacobiTheta2()
    t.testJacobiTheta3()
    t.testJacobiTheta4()
    t.test_JacobiTheta_complex()
    t.testJacobiEllipticSn()
    t.testJacobiEllipticCn()
    t.testJacobiEllipticDn()
    t.test_JacobiElliptic_complex()

if __name__ == '__main__':          # if run as script, run tests
    #test_elliptic_functions()
    unittest.main()



