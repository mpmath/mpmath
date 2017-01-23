class InverseLaplaceTransform(object):
    r"""
    Inverse Laplace transform methods are implemented using this
    class, in order to simplify the code and provide a common
    infrastructure.

    Implement a custom inverse Laplace transform algorithm by
    subclassing :class:`InverseLaplaceTransform` and implementing the
    appropriate methods. The subclass can then be used by
    :func:`~mpmath.invertlaplace` by passing it as the *method*
    argument.
    """

    def __init__(self,ctx):
        self.ctx = ctx

    def calc_laplace_parameter(self,t,**kwargs):
        """
        Determine the vector of Laplace parameter values needed for an
        algorithm, this will depend on the choice of algorithm (de
        Hoog is default), the algorithm-specific parameters passed (or
        default ones), and desired time.
        """
        raise NotImplementedError

    def calc_time_domain_solution(self,fp):
        """
        Compute the time domain solution, after computing the
        Laplace-space function evaluations at the abscissa required
        for the algorithm. Abscissa computed for one algorithm are
        typically not useful for another algorithm.
        """
        raise NotImplementedError

class FixedTalbot(InverseLaplaceTransform):

    def calc_laplace_parameter(self,t,**kwargs):
        r"""
        The "fixed" Talbot method deforms the Bromwich contour towards
        `-\infty` in the shape of a parabola. Traditionally the Talbot
        algorithm has adjustable parameters, but the "fixed" version
        does not. The `r` parameter could be passed in as a parameter,
        if you want to override the default given by (Abate & Valko,
        2004).

        The Laplace parameter is sampled along a parabola opening
        along the negative imaginary axis, with the base of the
        parabola along the real axis at
        `p=\frac{r}{t_\mathrm{max}}`. As the number of terms used in
        the approximation (degree) grows, the abscissa required for
        function evaluation tend towards `-\infty`, requiring high
        precision to prevent overflow.

        **Optional arguments**

        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        recognizes the following keywords

        *tmax*
            maximum time associated with vector of times
            (typically just the time requested)
        *degree*
            integer order of approximation (M = number of terms)
        *r*
            abscissa for `p_0` (otherwise computed using rule
            of thumb `2M/5`)

        The working precision will be increased according to a rule of
        thumb. If 'degree' is not specified, the working precision and
        degree are chosen to hopefully achieve the dps of the calling
        context. If 'degree' is specified, the working precision is
        chosen to achieve maximum resulting precision for the
        specified degree.

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
        # maximum time desired (used for scaling) default is requested
        # time.
        self.tmax = self.ctx.convert(kwargs.get('tmax',self.t))

        # empirical relationships used here based on a linear fit of
        # requested and delivered dps for exponentially decaying time
        # functions for requested dps up to 512.

        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = self.degree
        else:
            self.dps_goal = int(1.72*self.ctx.dps)
            self.degree = max(12,int(1.38*self.dps_goal))

        M = self.degree

        # this is adjusting the dps of the calling context hopefully
        # the caller doesn't monkey around with it between calling
        # this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = self.dps_goal

        # Abate & Valko rule of thumb for r parameter
        self.r = kwargs.get('r',self.ctx.fraction(2,5)*M)

        self.theta = self.ctx.linspace(0.0, self.ctx.pi, M+1)

        self.cot_theta = self.ctx.matrix(M,1)
        self.cot_theta[0] = 0 # not used

        # all but time-dependent part of p
        self.delta = self.ctx.matrix(M,1)
        self.delta[0] = self.r

        for i in range(1,M):
            self.cot_theta[i] = self.ctx.cot(self.theta[i])
            self.delta[i] = self.r*self.theta[i]*(self.cot_theta[i] + 1j)

        self.p = self.ctx.matrix(M,1)
        self.p = self.delta/self.tmax

        # NB: p is complex (mpc)

    def calc_time_domain_solution(self,fp,t,manual_prec=False):
        r"""
        The fixed Talbot time-domain solution is computed from the
        Laplace-space function evaluations using

        `f(t,M)=\frac{2}{5t}\sum_{k=0}^{M-1}\mathrm{Re} \left[
        \gamma_k \bar{f}(p_k)\right]`

        where

        `\gamma_0 = \frac{1}{2}e^{r}\bar{f}(p_0)`

        `\gamma_k = e^{tp_k}\left\lbrace 1 + \frac{jk\pi}{M}\left[1 +
        \cot \left( \frac{k \pi}{M} \right)^2 \right] - j\cot\left(
        \frac{k \pi}{M}\right)\right \rbrace \qquad 1\le k<M`.

        Again, `j=\sqrt{-1}`.

        Before calling this function, call
        :class:`~mpmath.calculus.inverselaplace.FixedTalbot.calc_laplace_parameter`
        to set the parameters and compute the required coefficients.
        """

        # required
        # ------------------------------
        self.t = self.ctx.convert(t)

        # assume fp was computed from p matrix returned from
        # calc_laplace_parameter(), so is already a list or matrix of
        # mpmath 'mpc' types

        # these were computed in previous call to
        # calc_laplace_parameter()
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

        # setting dps back to value when calc_laplace_parameter was
        # called, unless flag is set.
        if not manual_prec:
            self.ctx.dps = self.dps_orig

        return result.real

# ****************************************

class Stehfest(InverseLaplaceTransform):

    def calc_laplace_parameter(self,t,**kwargs):
        r"""
        The Gaver-Stehfest method is a discrete approximation of the
        Widder-Post inversion algorithm, rather than a direct
        approximation of the Bromwich contour integral.

        The method abscissa along the real axis, and therefore has
        issues inverting oscillatory functions (which have poles in
        pairs away from the real axis).

        The working precision will be increased according to a rule of
        thumb. If 'degree' is not specified, the working precision and
        degree are chosen to hopefully achieve the dps of the calling
        context. If 'degree' is specified, the working precision is
        chosen to achieve maximum resulting precision for the
        specified degree.

        `p_j = \frac{j \log 2}{t} \qquad 1 \le j \le M`
        """

        # required
        # ------------------------------
        # time of desired approximation
        self.t = self.ctx.convert(t)

        # optional
        # ------------------------------

        # empirical relationships used here based on a linear fit of
        # requested and delivered dps for exponentially decaying time
        # functions for requested dps up to 512.

        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = int(1.38*self.degree)
        else:
            self.dps_goal = int(2.93*self.ctx.dps)
            self.degree = max(16,self.dps_goal)

        # _coeff routine requires even degree
        if self.degree%2 > 0:
            self.degree += 1

        M = self.degree

        # this is adjusting the dps of the calling context
        # hopefully the caller doesn't monkey around with it
        # between calling this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = self.dps_goal

        self.V = self._coeff()
        self.p = self.ctx.matrix(self.ctx.arange(1,M+1))*self.ctx.ln2/self.t

        # NB: p is real (mpf)

    def _coeff(self):
        r"""Salzer summation weights (aka, "Stehfest coefficients")
        only depend on the approximation order (M) and the precision"""

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

    def calc_time_domain_solution(self,fp,t,manual_prec=False):
        r"""
        Compute time-domain Stehfest algorithm solution.

        `f(t,M) = \frac{\log 2}{t} \sum_{k=1}^{M} V_k \bar{f}\left(
        p_k \right)`

        where

        `V_k=(-1)^{k + N/2}\sum^{\min(k,N/2)}_{j=\lfloor(k+1)/2
        \rfloor}\frac{j^{\frac{N}{2}}(2j)!}{(\frac{N}{2}-j)!\,j!\,
        (j-1)!\,(k-j)!\,(2j-k)!}.
        """

        # required
        self.t = self.ctx.convert(t)

        # assume fp was computed from p matrix returned from
        # calc_laplace_parameter(), so is already
        # a list or matrix of mpmath 'mpf' types

        result = self.ctx.fdot(self.V,fp)*self.ctx.ln2/self.t

        # setting dps back to value when calc_laplace_parameter was called
        if not manual_prec:
            self.ctx.dps = self.dps_orig

        # ignore any small imaginary part
        return result.real

# ****************************************

class deHoog(InverseLaplaceTransform):

    def calc_laplace_parameter(self,t,**kwargs):
        r"""the de Hoog, Knight & Stokes algorithm is an
        accelerated form of the Fourier series numerical
        inverse Laplace transform algorithms, like Crump (YEAR)
        or Dubner & Abate (YEAR).

        `p = `
        """

        self.t = self.ctx.convert(t)
        self.tmax = self.t

        # empirical relationships used here based on a linear fit of
        # requested and delivered dps for exponentially decaying time
        # functions for requested dps up to 512.

        if 'degree' in kwargs:
            self.degree = kwargs['degree']
            self.dps_goal = int(1.38*self.degree)
        else:
            self.dps_goal = int(self.ctx.dps*1.36)
            self.degree = max(10,self.dps_goal)

        # 2*M+1 terms in approximation
        M = self.degree

        # adjust alpha component of abscissa of convergence for higher
        # precision
        tmp = self.ctx.power(10.0,-self.dps_goal)
        self.alpha = self.ctx.convert(kwargs.get('alpha',tmp))

        # desired tolerance (here simply related to alpha)
        self.tol = self.ctx.convert(kwargs.get('tol',self.alpha*10.0))
        self.np = 2*self.degree+1 # number of terms in approximation

        # this is adjusting the dps of the calling context
        # hopefully the caller doesn't monkey around with it
        # between calling this routine and calc_time_domain_solution()
        self.dps_orig = self.ctx.dps
        self.ctx.dps = self.dps_goal

        # scaling factor (likely tun-able, but 2 is typical)
        self.scale = kwargs.get('scale',2)
        self.T = self.ctx.convert(kwargs.get('T',self.scale*self.tmax))

        self.p = self.ctx.matrix(2*M+1,1)
        self.gamma = self.alpha - self.ctx.log(self.tol)/(self.scale*self.T)
        self.p = (self.gamma + self.ctx.pi*
                  self.ctx.matrix(self.ctx.arange(self.np))/self.T*1j)

        # NB: p is complex (mpc)

    def calc_time_domain_solution(self,fp,t,manual_prec=False):
        r"""Calculate time-domain solution for
        de Hoog, Knight & Stokes algorithm.

        The un-accelerated Fourier series approach:

        `f(t,2M+1) = \frac{e^{\gamma t}}{T} \sum_{k=0}^{2M}{}^{'}
        \Re\left[\bar{f}\left(\gamma_0 +\frac{i\pi k}{T}\right)
        \exp\left(\frac{i\pi t}{T}\right)\right],` 
        """
        
        M = self.degree
        np = self.np
        T = self.T

        self.t = self.ctx.convert(t)

        # would it be useful to try re-using
        # space between e&q and A&B?
        e = self.ctx.zeros(np,M+1)
        q = self.ctx.matrix(np,M)
        d = self.ctx.matrix(np,1)
        A = self.ctx.zeros(np+2,1)
        B = self.ctx.ones(np+2,1)

        # initialize Q-D table
        # e[0:2*M,0] = 0.0 + 0.0j
        q[0,0] = fp[1]/(fp[0]/2)
        for i in range(1,2*M):
            q[i,0] = fp[i+1]/fp[i]

        # rhombus rule for filling triangular Q-D table (e & q)
        for r in range(1,M+1):
            # start with e, column 1, 0:2*M-2
            mr = 2*(M-r)
            e[0:mr,r] = q[1:mr+1,r-1] - q[0:mr,r-1] + e[1:mr+1,r-1]
            if not r == M:
                rq = r+1
                mr = 2*(M-rq)+1
                for i in range(mr):
                    q[i,rq-1] = q[i+1,rq-2]*e[i+1,rq-1]/e[i,rq-1]

        # build up continued fraction coefficients (d)
        d[0] = fp[0]/2
        for r in range(1,M+1):
            d[2*r-1] = -q[0,r-1] # even terms
            d[2*r]   = -e[0,r]   # odd terms

        # seed A and B for recurrence
        #A[0] = 0.0 + 0.0j
        A[1] = d[0]
        #B[0:2] = 1.0 + 0.0j

        # base of the power series
        z = self.ctx.expjpi(self.t/T) # i*pi is already in fcn

        # coefficients of Pade approximation (A & B)
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
        result = self.ctx.exp(self.gamma*self.t)/T*(A[np]/B[np]).real

        # setting dps back to value when calc_laplace_parameter was called
        if not manual_prec:
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
        evaluated is assumed to be a real-valued function of time.

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
        >>> ft = lambda t: t*exp(-t)
        >>> [(t,ft(t),ft(t)-invertlaplace(fp,t,method='talbot')) for t in tvec]
        [(0.001, mpf('0.0009990004998333751'), mpf('8.5792304356121216e-20')),
        (0.01, mpf('0.0099004983374916811'), mpf('3.2700764669804695e-19')),
        (0.1, mpf('0.09048374180359596'), mpf('-1.7521580005216785e-18')),
        (1.0, mpf('0.36787944117144233'), mpf('1.2428864009344032e-17')),
        (10.0, mpf('0.00045399929762484856'), mpf('4.0451348930665763e-20'))]

        The methods also work for higher precision:

        >>> mp.dps = 100
        >>> [(t,nstr(ft(t),15),nstr(ft(t)-invertlaplace(fp,t,method='talbot'),15)) for t in tvec]
        [(0.001, '0.000999000499833375', '-4.96868310693356e-105'),
        (0.01, '0.00990049833749168', '1.23032291513122e-104'),
        (0.1, '0.090483741803596', '7.25752777054208e-103'),
        (1.0, '0.367879441171442', '1.0212982518952e-102'),
        (10.0, '0.000453999297624849', '1.43644942538356e-105')]

        `\bar{f}(p) = \frac{1}{p^2+1}`

        `f(t) = \mathrm{J}_0(t)`

        >>> mp.dps = 15
        >>> fp = lambda p: 1/sqrt(p*p + 1)
        >>> ft = lambda t: besselj(0,t)
        >>> [(t,ft(t),ft(t)-invertlaplace(fp,t)) for t in tvec]
        [(0.001, mpf('0.99999975000001562'), mpf('-8.2477943034013953e-18')),
        (0.01, mpf('0.99997500015624957'), mpf('-3.6981014489887214e-17')),
        (0.1, mpf('0.99750156206604003'), mpf('1.7277947790359516e-18')),
        (1.0, mpf('0.76519768655796655'), mpf('1.0619604044683858e-17')),
        (10.0, mpf('-0.24593576445134834'), mpf('-3.2166963159814471e-17'))]

        `\bar{f}(p) = \frac{\log p}{p}`

        `f(t) = -\gamma -\log t`

        >>> fp = lambda p: log(p)/p
        >>> ft = lambda t: -euler-log(t)
        >>> [(t,ft(t),ft(t)-invertlaplace(fp,t,method='stehfest')) for t in tvec]
        [(0.001, mpf('6.3305396140806041'), mpf('-1.9212663483786344e-16')),
        (0.01, mpf('4.0279545210865582'), mpf('-4.814860932007035e-16')),
        (0.1, mpf('1.7253694280925127'), mpf('-7.2309776669517614e-18')),
        (1.0, mpf('-0.57721566490153287'), mpf('2.4823521988792137e-17')),
        (10.0, mpf('-2.8798007578955787'), mpf('-1.1353030106517025e-16'))]

        **Options**

        :func:`~mpmath.invertlaplace` recognizes the following keywords
        valid for all methods:

        *method*
            Chooses numerical inverse Laplace transform algorithm
            (described below).
        *degree*
            Number of terms used in the approximation

        **Algorithms**

        Mpmath implements three numerical inverse Laplace transform
        algorithms, attributed to: Talbot, Stehfest, and de Hoog,
        Knight and Stokes. These can be selected by using
        *method='talbot'*, *method='stehfest'*, or *method='dehoog'*
        or by passing the classes *method=FixedTalbot*,
        *method=Stehfest*, or *method=deHoog*. The functions
        :func:`~mpmath.invlaptalbot`, :func:`~mpmath.invlapstehfest`,
        and :func:`~mpmath.invlapdehoog` are also available as
        shortcuts.

        All three algorithms implement a heuristic balance between the
        requested precision, and the precision used internally for the
        calculations. This has been tuned for precision up to few
        hundred decimal digits. The Laplace transform converts
        behavior for time (i.e., along a line) into the entire complex
        `p`-plane.  Singularities, poles, and branch cuts in the
        complex `p`-plane contain all the information regarding the
        time behavior of the function of interest. Any numerical
        method must therefore sample `p`-plane "close enough" to the
        singularities to accurately characterize them, while not
        getting too close to have catastrophic cancellation issues. If
        one or more of the singularities in the `p`-plane is not on
        the left side of the Bromwich contour, its effects will be
        left out of the computed solution, and the answer will be
        wrong.

        The fixed Talbot method is high accuracy and fast, but the
        method can catastrophically fail for certain time-domain
        behaviors, including a Heaviside step function for positive
        time (e.g., `H(t-2)`), or some oscillatory behaviors. The
        Talbot method usually has adjustable parameters, but the
        "fixed" variety implemented here does not.  This function is
        not defined for small time (in this case `t<2`).  This method
        deforms the Bromwich integral contour in the shape of a
        parabola headed towards `-\infty`, which leads to problems
        when the solution has a decaying exponential in it (e.g., a
        Heaviside step function is equivalent to multiplying by a
        decaying exponential in Laplace space).

        The Stehfest algorithm only uses abscissa along the real axis
        of the imaginary plane to estimate the time-domain
        function. Oscillatory time-domain functions have poles away
        from the real axis, so this method does not work well with
        oscillatory functions. This method also depends on summation
        of terms in a series that grow very large, and will have
        catastrophic cancellation during summation if the working
        precision is too low.

        The de Hoog, Knight and Stokes method is essentially a
        Fourier-series quadrature type approximation to the Bromwich
        contour integral, with non-linear series acceleration and an
        analytical expression for the remainder term. This method is
        typically the most robust and is therefore the default method.

        **Singularities**

        All numerical inverse Laplace transform methods have problems
        at large time when the Laplace-space function has poles,
        singularities, or branch cuts to the right of the origin in
        the complex plane. For simple poles in `\bar{f}(p)`, at the
        `p`-plane origin is a constant in time. A pole to the left of
        the origin is a decreasing function of time, and a pole to the
        right of the origin leads to an increasing function in time
        (see Duffy [3] Figure 4.10.4, p. 228).

        For example:

        `\bar{f}(p)=\frac{1}{p^2-9}`

        `f(t)=\frac{1}{3}\sinh 3t`

        In general as `p \rightarrow \infty` `t \rightarrow 0` and
        vice-versa. All numerical inverse Laplace transform methods
        require their abscissa to shift closer to the origin for larger
        times. If the abscissa shift left of the rightmost singularity
        in the Laplace domain, the answer will be completely wrong.

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
        # at the required abscissa.
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
