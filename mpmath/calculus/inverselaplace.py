
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
        self.degree = None
        
        # number of digits to add to working precision
        # really just a guess (could be tuned for each method)
        self.step = max(20,int(0.2*self.ctx.dps))

        # cache some coefficients for talbot and stehfest methods
        self.talbot_cache = {}
        self.stehfest_cache = {}
        # no caching for deHoog et al. method
        
    def clear(self):
        """
        Delete cached coefficient data.
        """
        self.talbot_cache = {}
        self.stehfest_cache = {}
        
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

        The laplace parameter is sampled along a parabola 
        opening along the negative imaginary axis, with the 
        bottom of the parabola at `p=\frac{r}{t_\mathrm{max}}+0j`

        **Optional arguments**
        
        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        recognizes the following keywords

        *tmax*
            maximum time associated with vector of times 
            (typically just the time requested)
        *degree*
            integer order of approximation (M = number of terms)
        *dps*
            desired decimal precision (otherwise computed from
            requested degree via heuristic). 
        *r*
            abcissa for `p_0` (otherwise computed using rule
            of thumb `2M/5`)
        
        The working precision will be increased according to a
        rule of thumb, in relation to the degree; dps_goal`=1.7M`. 
        dps_goal can be overridden with a user-specified value.

        `p_0=\frac{r}{t}`

        `p_i=\frac{i r \pi}{Mt_\mathrm{max}}\left[\cot\left(
        \frac{i\pi}{M}\right) + j \right] \qquad 1\le i <M`

        where `j=\sqrt{-1}`, `r=2M/5`, and `t_\mathrm{max}` is the 
        maximum specified time.
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

        # order of approximation
        # rule for extended precision from Abate & Valko (2004)
        # "Multi-precision Laplace Transform Inversion"
        
        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = max(self.ctx.dps, int(1.7*self.degree))
        else:
            self.degree = int(self.ctx.dps/1.7)
            self.dps_goal = self.ctx.dps

        self.step = kwargs.get('step',self.step)
            
        M = self.degree
        dps_inner = self.dps_goal + self.step

        # this is adjusting the dps of the calling context
        # hopefully the caller doesn't monkey around with it
        # between calling this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = dps_inner
        
        # Abate & Valko rule of thumb for r parameter
        self.r = kwargs.get('r',self.ctx.fraction(2,5)*M)
        
        if (dps_inner,M) in self.talbot_cache:
            (self.theta, self.cot_theta, self.delta) = self.talbot_cache[dps_inner,M]

        else:
            self.theta = self.ctx.linspace(0.0, self.ctx.pi, M+1)
                            
            self.cot_theta = self.ctx.matrix(M,1) 
            self.cot_theta[0] = 0 # not used

            # all but time-dependent part of p
            self.delta = self.ctx.matrix(M,1) 
            self.delta[0] = self.r
            
            for i in range(1,M):
                self.cot_theta[i] = self.ctx.cot(self.theta[i])
                self.delta[i] = self.r*self.theta[i]*(self.cot_theta[i] + 1j)
                
            self.talbot_cache[dps_inner,M] = (self.theta, self.cot_theta, self.delta)

        # don't cache p; it depends on t_max
        self.p = self.ctx.matrix(M,1)

        for i in range(M):
            self.p[i] = self.delta[i]/self.tmax

        # NB: p is complex (mpc)
    
    def calc_time_domain_solution(self,fp,t):
        r"""
        The time-domain solution is computed from the Laplace-space
        function evaluations using

        `f(t,M)=\frac{2}{5t}\sum_{k=0}^{M-1}\mathrm{Re} 
        \left[ \gamma_k \bar{f}(p_k)\right]`

        where 

        `\gamma_0 = \frac{1}{2}e^{r}\bar{f}(p_0)`

        `\gamma_k = e^{tp_k}\left\lbrace 1 + \frac{jk\pi}{M}\left[1 + 
        \cot \left( \frac{k \pi}{M} \right)^2 \right] - 
        j\cot\left( \frac{k \pi}{M}\right)\right \rbrace \qquad 1\le k<M`.

        Again, `j=\sqrt{-1}`.

        Before calling this function, call 
        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        to set the parameters and compute the required coefficients.
        """

        # required
        # ------------------------------
        self.t = self.ctx.convert(t)
        
        # assume fp was computed from p matrix returned from
        # calc_laplace_parameter(), so is already
        # a list or matrix of mpmath 'mpc' types

        # these were computed in previous call to calc_laplace_parameter()
        theta = self.theta
        delta = self.delta
        M = self.degree
        p = self.p
        r = self.r
        
        ans = self.ctx.matrix(M,1)
        ans[0] = self.ctx.exp(delta[0])*fp[0]/2
        
        for i in range(1,M):
            ans[i] = self.ctx.exp(delta[i])*fp[i]*(
                1 + 1j*theta[i]*(1 + self.cot_theta[i]**2) -
                1j*self.cot_theta[i])

        result = self.ctx.fraction(2,5)*self.ctx.fsum(ans)/self.t

        # setting dps back to value when calc_laplace_parameter was called
        self.ctx.dps = self.dps_orig
        
        return result.real

# ****************************************

class Stehfest(InverseLaplaceTransform):

    def calc_laplace_parameter(self, t, **kwargs):
        r"""
        The Gaver-Stehfest method is a discrete approximation of the
        Widder-Post inversion algorithm, rather than a direct
        approximation of the Bromwich contour integral.

        The method only uses values of `p` along the real axis, and
        therefore has issues inverting oscillatory functions (which
        have poles in pairs away from the real axis).

        The Stehfest coefficients only depend on the order of
        approximation and not the desired time, so they could be
        computed once and re-used for a given degree and working
        precision.
        """

        # required
        # ------------------------------
        # time of desired approximation
        self.t = self.ctx.convert(t)

        # optional
        # ------------------------------
        
        # rule-of-thumb for extended precision from Abate & Valko (2004)
        # "Multi-precision Laplace Transform Inversion"

        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = max(self.ctx.dps,int(2.1*self.degree))
        else:
            self.degree = int(self.ctx.dps/2.1)
            self.dps_goal = self.ctx.dps
                    
        # _coeff routine requires degree must be even
        if self.degree%2 > 0:
            self.degree += 1

        self.step = kwargs.get('step',self.step)
            
        M = self.degree
        dps_inner = self.dps_goal + self.step

        # this is adjusting the dps of the calling context
        # hopefully the caller doesn't monkey around with it
        # between calling this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = dps_inner
        
        if (dps_inner,M) in self.stehfest_cache:
            self.V = self.stehfest_cache[dps_inner,M]
            
        else:
            self.V = self._coeff()
            self.stehfest_cache[dps_inner,M] = self.V
                
        # p not cached since it depends on t, and is trivial to compute
        self.p = self.ctx.matrix(self.ctx.arange(1,M+1))*self.ctx.ln2/self.t

        # NB: p is real (mpf)
        
    def _coeff(self):
        """The Stehfest coefficients only depend on the 
        approximation order, M and precsion"""

        M = self.degree
        M2 = M/2

        V = self.ctx.matrix(M,1)

        # Salzer summation weights
        # get very large in magnitude and oscillate in sign,
        # if the precision is not high enough, there will be
        # catastrophic cancellation
        for k in range(1,M+1):
            z = self.ctx.matrix(int(min(k,M2)+1),1)
            for j in range(int(self.ctx.floor((k+1)/2.0)),min(k,M2)+1):
                z[j] = (self.ctx.power(j,M2)*self.ctx.fac(2*j)/
                        (self.ctx.fac(M2-j)*self.ctx.fac(j)*
                         self.ctx.fac(j-1)*self.ctx.fac(k-j)*
                         self.ctx.fac(2*j-k)))
            V[k-1] = self.ctx.power(-1,k+M2)*self.ctx.fsum(z)

        return V

    def calc_time_domain_solution(self,fp,t):
        """Compute time-domain solution using f(p)
        and coefficients"""

        # required
        self.t = self.ctx.convert(t)
        
        # assume fp was computed from p matrix returned from
        # calc_laplace_parameter(), so is already
        # a list or matrix of mpmath 'mpf' types

        result = self.ctx.fdot(self.V,fp)*self.ctx.ln2/t

        # setting dps back to value when calc_laplace_parameter was called
        self.ctx.dps = self.dps_orig
        
        # ignore any small imaginary part
        return result.real

# ****************************************

class deHoog(InverseLaplaceTransform):

    def calc_laplace_parameter(self, t, **kwargs):

        self.t = self.ctx.convert(t)

        # 2*M+1 terms used below

        # no widely-used rule of thumb for increasing precsion as 
        # degree of approximation increases (based on experimentation)
        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = max(self.ctx.dps,int(1.8*self.degree))
        else:
            self.degree = int(self.ctx.dps/1.8)
            self.dps_goal = self.ctx.dps
            
        self.step = kwargs.get('step',self.step)
            
        M = self.degree
        dps_inner = self.dps_goal + self.step

        # heuristic based on desired precision
        tmp = self.ctx.power(10.0,-self.dps_goal)
        self.alpha = self.ctx.convert(kwargs.get('alpha',tmp))

        # desired tolerance
        self.tol = self.ctx.convert(kwargs.get('tol',self.alpha*10.0))
        self.np = 2*self.degree+1
        
        # this is adjusting the dps of the calling context
        # hopefully the caller doesn't monkey around with it
        # between calling this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = dps_inner
        
        # scaling factor (likely tunable, but 2 is typical)
        self.scale = kwargs.get('scale',2)
        self.T = self.ctx.convert(kwargs.get('T',self.scale*self.t))        

        # not chached since p depends on, alpha, tol, T, and scale,
        # and is not difficult to compute.
            
        self.p = self.ctx.matrix(2*M+1,1)
        self.gamma = self.alpha - self.ctx.log(self.tol)/(self.scale*self.T)
        for i in range(2*M+1):
            self.p[i] = self.gamma + self.ctx.pi*i/self.T*1j

        # NB: p is complex (mpc)

    def calc_time_domain_solution(self,fp,t):

        M = self.degree
        np = self.np
        T = self.T

        self.t = self.ctx.convert(t)

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
        result = self.ctx.exp(self.gamma*t)/T*(A[np]/B[np]).real

        # setting dps back to value when calc_laplace_parameter was called
        self.ctx.dps = self.dps_orig
        
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

        `\mathcal{L}\left[ f(t) \right]=\bar{f}(p)`

        `\mathcal{L}^{-1}\left[ \bar{f}(p) \right] = f(t)`

        `\bar{f}(p) = \frac{1}{(p+1)^2}`

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
        a parabola headed towards `-\infty`, which leads to problems when 
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

        1. NIST Digital Library of Mathematical Functions. Table 1.14.4 (http://dlmf.nist.gov/1.14T4)
        2. Cohen, A.M. (2007). Numerical Methods for Laplace Transform Inversion, Springer.
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
                raise ValueError("unknown invlap algorithm: %s" % rule)
        else:
            rule = rule(ctx)

        # determine the vector of Laplace-space parameter
        # needed for the requested method and desired time
        rule.calc_laplace_parameter(t,**kwargs)

        # compute the Laplace-space function evalutations
        # at the required abcissa.
        fp = map(f,rule.p)

        # compute the time-domain solution from the
        # Laplace-space function evaluations
        return rule.calc_time_domain_solution(fp,t)

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
