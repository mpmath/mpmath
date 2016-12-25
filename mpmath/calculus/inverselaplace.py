
class InverseLaplaceTransform(object):
    r"""
    Inverse Lapalce Transform method are implemented using this class, 
    in order to simplify the code and provide a common infrastructure.

    You can implement a custom inverse Laplace transform algorithm by subclassing
    :class:`InverseLaplaceTransform` and implementing the appropriate
    methods. The subclass can then be used by :func:`~mpmath.invertlaplace` by
    passing it as the *method* argument.
    """

    def __init__(self, ctx):
        self.ctx = ctx

        # a default level of approximation appropriate  
        # for the given precision goal
        self.degree = 18

        # decimal digits of precision for computing p vectors, 
        # and f(t) solutions.  (Talbot and Stehfest have 
        # their own expressions that override this minimum)
        self.dps_goal = 30
        
    def calc_laplace_parameter(self,t,**kwargs):
        """
        Determine the vector of Laplace parameter values
        needed for an algorithm, this will depend on the choice
        of algorithm (de Hoog is default), the algorithm-specific
        parameters passed (or default ones), and desired time.
        """
        raise NotImplementedError

    def calc_time_domain_solution(self,fp):
        """
        Compute the time domain solution, after computing the
        Laplace-space function evalutations at the abcissa 
        required for the algorithm. Abcissa computed for one
        algorithm are typically not useful for another algorithm.
        """
        raise NotImplementedError

class FixedTalbot(InverseLaplaceTransform):
    
    def calc_laplace_parameter(self, t, **kwargs):
        r"""
        The "fixed" Talbot method deforms the Bromwich contour
        towards `-\infty` in the shape of a parabola. Traditionally
        the Talbot algorithm has adjustable parameters, but the
        "fixed" version does not. The `r` parameter could be 
        passed in as a parameter, if you want to override
        the default given by (Abate & Valko, 2004).
        """

        # required
        # ------------------------------
        # time of desired approximation
        self.t = self.ctx.convert(t)

        # optional
        # ------------------------------
        # maximum time desired (used for scaling)
        # default is requested time. 
        self.tmax = self.ctx.convert(kwargs.get('tmax',self.t))

        # integer order of approximation
        self.degree = kwargs.get('degree',self.degree)
        M = self.degree

        # rule for extended precision from Abate & Valko (2004)
        # "Multi-precision Laplace Transform Inversion"
        self.dps_goal = kwargs.get('dps',max(self.dps_goal,0.6*M))

        # Abate & Valko rule of thumb for r parameter
        self.r = kwargs.get('r',self.ctx.fraction(2,5)*M)
        
        self.p = self.ctx.matrix(M,1)
        self.cot_theta = self.ctx.matrix(M,1) # needed several times
        self.p[0] = self.r/self.tmax
        self.cot_theta[0] = 0.0 # not used

        with self.ctx.workdps(self.dps_goal):
            self.theta = self.ctx.linspace(0.0, self.ctx.pi, M+1)
            
            for i in range(1,M):
                self.cot_theta[i] = self.ctx.cot(self.theta[i])
                self.p[i] = self.r*self.theta[i]*( 
                    self.cot_theta[i] + 1j)/self.tmax
    
    def calc_time_domain_solution(self,fp):

        theta = self.theta
        M = self.degree
        t = self.t
        p = self.p
        r = self.r
        ans = self.ctx.matrix(M,1)

        with self.ctx.workdps(self.dps_goal):

            ans[0] = self.ctx.exp(r)*fp[0]/2
            
            for i in range(1,M):
                ans[i] = self.ctx.exp(t*p[i])*fp[i]*(
                    1 + 1j*theta[i]*(1 + self.cot_theta[i]**2) -
                    1j*self.cot_theta[i])
    
            result = self.ctx.fraction(2,5)*self.ctx.fsum(ans)/self.t

        return result.real

# ****************************************

class Stehfest(InverseLaplaceTransform):

    def calc_laplace_parameter(self, t, **kwargs):

        # time of desired approximation
        self.t = self.ctx.convert(t)

        self.degree = kwargs.get('degree',self.degree)

        # 2.1 = rule for extended precision from Abate & Valko (2004)
        # "Multi-precision Laplace Transform Inversion"
        self.dps_goal = kwargs.get('dps',max(self.dps_goal,2.1*self.degree))

        M = self.degree

        # don't compute V here
        self.V = kwargs.get('V',None)       

        with self.ctx.workdps(self.dps_goal):

            self.p = (self.ctx.matrix(self.ctx.arange(1,M+1))*
                      self.ctx.ln2/self.t)
    
    def _coeff(self):
        """Stehfest coefficients only depend on M"""

        M = self.degree
        # order must be odd
        if M%2 > 0:
            self.degree += 1
            M = self.degree
            
        M2 = M/2

        self.V = self.ctx.matrix(M,1)

        # Salzer summation weights
        for k in range(1,M+1):
            z = self.ctx.matrix(int(min(k,M2)+1),1)
            for j in range(int(self.ctx.floor((k+1)/2.0)),min(k,M2)+1):
                z[j] = (self.ctx.power(j,M2)*self.ctx.fac(2*j)/
                        (self.ctx.fac(M2-j)*self.ctx.fac(j)*
                         self.ctx.fac(j-1)*self.ctx.fac(k-j)*self.ctx.fac(2*j-k)))
            self.V[k-1] = self.ctx.power(-1,k+M2)*self.ctx.fsum(z)

    def calc_time_domain_solution(self,fp):
        """Compute time-domain solution using f(p)
        and coefficients"""

        with self.ctx.workdps(self.dps_goal):

            if self.V is None:
                self._coeff()
            else:
                self.V = self.ctx.convert(self.V)
    
            result = self.ctx.fdot(self.V,fp)*self.ctx.ln2/self.t

        # ignore any small imaginary part
        return result.real

# ****************************************

class deHoog(InverseLaplaceTransform):

    def calc_laplace_parameter(self, t, **kwargs):

        self.t = self.ctx.convert(t)

        # 2*M+1 terms used below
        self.degree = kwargs.get('degree',self.degree)

        # heuristic based on desired precision
        tmp = self.ctx.power(10.0,-self.dps_goal)
        self.alpha = self.ctx.convert(kwargs.get('alpha',tmp))

        # desired tolerance
        self.tol = self.ctx.convert(kwargs.get('tol',self.alpha*10.0))
        self.np = 2*self.degree+1

        # scaling factor (this is likely tunable)
        self.scale = kwargs.get('scale',2)
        self.T = self.ctx.convert(kwargs.get('T',self.scale*self.t))        

        # no widely-used rule of thumb for increasing precsion as 
        # degree of approximation increases, this based on limited experience
        self.dps_goal = kwargs.get('dps',max(self.dps_goal,1.8*self.degree))

        M = self.degree
        T = self.T
        self.p = self.ctx.matrix(2*M+1,1)

        with self.ctx.workdps(self.dps_goal):

            for i in range(2*M+1):
                self.p[i] = (self.alpha - 
                             self.ctx.log(self.tol)/(self.scale*T) + 
                             self.ctx.pi*i/T*1j)

    def calc_time_domain_solution(self,fp):

        M = self.degree
        np = self.np
        t = self.t
        T = self.T
        alpha = self.alpha
        tol = self.tol

        with self.ctx.workdps(self.dps_goal):
            
            gamma = alpha - self.ctx.log(tol)/(self.scale*T)

            e = self.ctx.matrix(np,M+1)
            q = self.ctx.matrix(np,M)
            d = self.ctx.matrix(np,1)
            A = self.ctx.matrix(np+2,1)
            B = self.ctx.matrix(np+2,1)
        
            # initialize Q-D table
            e[0:2*M,0] = 0.0
            q[0,0] = fp[1]/(fp[0]/2)
            for i in range(1,2*M):
                q[i,0] = fp[i+1]/fp[i]

            # rhombus rule for filling triangular Q-D table
            for r in range(1,M+1):
                # start with e, column 1, 0:2*M-2
                mr = 2*(M-r)
                e[0:mr,r] = q[1:mr+1,r-1] - q[0:mr,r-1] + e[1:mr+1,r-1]
                if not r == M:
                    rq = r+1
                    mr = 2*(M-rq)+1
                    for i in range(mr):
                        q[i,rq-1] = q[i+1,rq-2]*e[i+1,rq-1]/e[i,rq-1]

            # build up continued fraction coefficients
            d[0] = fp[0]/2
            for r in range(1,M+1):
                d[2*r-1] = -q[0,r-1] # even terms
                d[2*r]   = -e[0,r]   # odd terms

            # seed A and B for recurrence
            A[0] = 0.0
            A[1] = d[0]
            B[0:2] = 1.0

            # base of the power series
            z = self.ctx.expjpi(t/T) # i*pi is already in fcn

            # coefficients of Pade approximation
            # using recurrence for all but last term
            for i in range(1,2*M):
                A[i+1] = A[i] + d[i]*A[i-1]*z
                B[i+1] = B[i] + d[i]*B[i-1]*z

            # "improved remainder" to continued fraction
            brem  = (1 + (d[2*M-1] - d[2*M])*z)/2
            # powm1(x,y) computes x^y - 1 more accurately near zero
            rem = brem*self.ctx.powm1(1 + d[2*M]*z/brem,
                                      self.ctx.fraction(1,2))

            # last term of recurrence using new remainder
            A[np] = A[2*M] + rem*A[2*M-1]
            B[np] = B[2*M] + rem*B[2*M-1]

            # diagonal Pade approximation
            # F=A/B represents accelerated trapezoid rule
            result = self.ctx.exp(gamma*t)/T*(A[np]/B[np]).real

        return result

# ****************************************

class LaplaceTransformInversionMethods:
    def __init__(ctx, *args, **kwargs):
        ctx._fixed_talbot = FixedTalbot(ctx)
        ctx._stehfest = Stehfest(ctx)
        ctx._de_hoog = deHoog(ctx)

    def invertlaplace(ctx, f, t, **kwargs):
        r"""
        Computes the numerical inverse Laplace transform for a 
        Laplace-space function at a given time.  The function being 
        evalutated is assumed to be a real-valued function of time.

        The user must supply a Laplace-space function (`\bar{f}(p)`), 
        and a desired time to estimate the time-domain solution (`f(t)`) 
        at.

        A few basic examples of Laplace-space functions with known 
        inverses (see references [1,2]) ::

        >>> from mpmath import *
        >>> mp.dps = 15
        >>> tvec = [0.001, 0.01, 0.1, 1, 10]       

        `\bar{f}(p) = \frac{1}{p+1}^2`
        `\mathcal{L}^{-1}\left\lbrace \bar{f}(p) \right\rbrace = f(t)`
        `f(t) = t e^{-t}`

        >>> fp = lambda p: 1/(p+1)**2
        >>> ft = labmda t: t*exp(-t)
        >>> [(t,ft(t),ft(t)-invertlaplace(fp,t)) for t in tvec]
        [(0.001, 0.000999000499833375, 3.66645177352691e-14),
        (0.01, 0.00990049833749168, 1.85540141845894e-14),
        (0.1, 0.090483741803596, 2.04666254170337e-12),
        (1.0, 0.367879441171442, 1.49428280490678e-11),
        (10.0, 0.000453999297624849, -1.34042555175552e-10)]

        `\bar{f}(p) = \frac{1}{p^2+1}`
        `f(t) = \mathrm{J}_0(t)`

        >>> fp = lambda p: 1/sqrt(p*p + 1)
        >>> ft = lambda t: besselj(0,t)
        >>> [(t,ft(t),ft(t)-invertlaplace(fp,t)) for t in tvec]
        [(0.001, 0.999999750000016, 7.65046554765067e-10),
        (0.01, 0.99997500015625, -7.68641303866616e-10),
        (0.1, 0.99750156206604, 2.5249371516145e-10),
        (1.0, 0.765197686557967, -2.35601046132965e-10),
        (10.0, -0.245935764451348, 8.00007627560062e-10)]

        `\bar{f}(p) = \frac{\log p}{p}`
        `f(t) = -\gamma -\log t`

        >>> fp = lambda p: log(p)/p
        >>> ft = lambda t: -euler-log(t)
        [(0.001, 6.3305396140806, 1.33936103728218e-9),
        (0.01, 4.02795452108656, 5.02676075666929e-10),
        (0.1, 1.72536942809251, -8.78157818978427e-10),
        (1.0, -0.577215664901533, 1.31568888918193e-9),
        (10.0, -2.87980075789558, 7.0671105477838e-10)]


        **Options**

        :func:`~mpmath.invertlaplace` recognizes the following keywords 
        valid for all methods:

        *method*
            Choses numerical inverse Laplace transform algorithm 
            (described below).
        *degree*
            Number of terms used in the approximation 
        *dps_goal*
            Goal for precision of the approximation

        **Algorithms**

        Mpmath implements three numerical inverse Laplace transform
        algorithms, attributed to: Talbot, Stehfest, and de Hoog, Knight 
        and Stokes. These can be selected by using *method='talbot'*, 
        *method='stehfest'*, or *method='dehoog'* or by passing the 
        classes *method=FixedTalbot*, *method=Stehfest*, or 
        *method=deHoog*. The functions :func:`~mpmath.invlaptalbot`, 
        :func:`~mpmath.invlapstehfest`, and :func:`~mpmath.invlapdehoog` 
        are also available as shortcuts.

        All three algorithms have the property that doubling the number of 
        evaluation points roughly doubles the accuracy (with some 
        additional limitations outlined below), so the methods work for 
        high precision. The Laplace transform converts behavior for time 
        (i.e., along a line) into the entire complex `p`-plane. 
        Singularities, poles, and branch cuts in the complex `p`-plane 
        contain all the information regarding the time behavior of the 
        function of interest. Any numerical method must therefore 
        sample `p`-plane "close enough" to the singularities to accurately 
        characterize them, while not getting too close to have catastrophic
        cancellation issues. If one or more of the singularities in the 
        `p`-plane is not on the left side of the Bromwich contour, its 
        effects will be left out of the computed solution, and the answer 
        will be wrong.

        The fixed Talbot method is high accuracy, except when the time-
        domain function has a Heaviside step function for positive time 
        (e.g., `H(t-2)`). The Talbot method usually has adjustable 
        parameters, but the "fixed" variety implemented here does not. 
        This function is not defined for small time (in this case `t<2`). 
        This method deforms the Bromwich integral contour in the shape of 
        a parabola headed towards `-infty`, which leads to problems when 
        the solution has a decaying exponential in it (e.g., a Heaviside 
        step function is equivalient to multiplying by a decaying
        exponential in Laplace space). 

        The Stehfest algorithm only uses abcissa along the real axis of the
        imaginary plane to estimate the time-domain function. Oscillatory 
        time-domain functions have poles away from the real axis, so this 
        method does not work well with oscillatory functions. This method 
        also depends on summation of terms in a series that grow very 
        large, and will have catastrophic cancellation during summation 
        if the  working precision is too low.

        The de Hoog, Knight and Stokes method is essentially a 
        Fourier-series quadratyre type approximation to the Bromwich 
        contour integral, with non-linear series acceleration and an 
        analytical expression for the remainder term. This method is 
        typically the most robust and is therefore the default method. 

        **Singularities**

        All numerical inverse Laplace transform methods have problems at 
        large time when the Laplace-space function has poles, 
        singularities, or branch cuts to the right of the origin in the 
        complex plane. For simple poles in `\bar{f}(p)`, at the `p`-plane 
        origin is a constant in time. A pole to the left of the origin is 
        a decreasing function of time, and a pole to the right of the 
        origin leads to an increasing function in time (see Duffy [3] 
        Figure 4.10.4, p. 228). 

        For example:

        `\bar{f}(p)=\frac{1}{p^2-9}`
        `f(t)=\frac{1}{3}\sinh 3t`

        In general as `p \rightarrow \infty` `t \rightarrow 0` and 
        vice-versa. All numerical inverse Laplace transform methods 
        require their abcissa to shift closer to the origin for larger 
        times. If the abcissa shift left of the rightmost singularity in 
        the Laplace domain, the answer will be completely wrong.

        **Lower-level Use**

        The algorithms are designed to also be used with numerical 
        methods, rather that just with simple analytical functions. The

        **References**

        1. NIST Digital Library of Mathematical Functions. Table 1.14.4 
             (http://dlmf.nist.gov/1.14T4)
        2. Cohen, A.M. (2007). Numerical Methods for Laplace Transform 
             Inversion, Springer.
        3. Duffy, D.G. (1998). Advanced Engineering Mathematics, CRC Press.

        """
        
        rule = kwargs.get('method','dehoog')
        if type(rule) is str:
            lrule = rule.lower()
            if lrule == 'talbot':
                rule = ctx._fixed_talbot
            elif lrule == 'stehfest':
                rule = ctx._stehfest
            elif lrule == 'dehoog':
                rule = ctx._de_hoog
            else:
                raise ValueError("unknown method: %s" % rule)
        else:
            rule = rule(ctx)

        # determine the vector of Laplace-space parameter
        # needed for the requested method and desired time
        rule.calc_laplace_parameter(t,**kwargs)

        # compute the Laplace-space function evalutations
        # at the required abcissa.
        
        np = rule.p.rows # p as column vector
        fp = ctx.matrix(np,1)
        for i in range(np):
            fp[i] = f(rule.p[i])

        # compute the time-domain solution from the
        # Laplace-space function evaluations
        v = rule.calc_time_domain_solution(fp)
        return v

    # shortcuts for the above function for specific methods 
    def invlaptalbot(ctx, *args, **kwargs):
        kwargs['method'] = 'talbot'
        return ctx.invertlaplace(*args, **kwargs)

    def invlapstehfest(ctx, *args, **kwargs):
        kwargs['method'] = 'stehfest'
        return ctx.invertlaplace(*args, **kwargs)

    def invlapdehoog(ctx, *args, **kwargs):
        kwargs['method'] = 'dehoog'
        return ctx.invertlaplace(*args, **kwargs)

# ****************************************

if __name__ == '__main__':
    import doctest
    doctest.testmod()
