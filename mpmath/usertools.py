
def monitor(f, input='print', output='print'):
    """
    Returns a wrapped copy of *f* that monitors evaluation by calling
    *input* with every input (*args*, *kwargs*) passed to *f* and
    *output* with every value returned from *f*. The default action
    (specify using the special string value ``'print'``) is to print
    inputs and outputs to stdout, along with the total evaluation
    count::

        >>> from mpmath import mp, diff, monitor, exp, findroot, sin
        >>> mp.dps = 5
        >>> diff(monitor(exp), 1)   # diff will eval f(x-h) and f(x+h)
        in  0 (mpf('0.9999999990686774253845215'),) {}
        out 0 mpf('2.718281825927448005528206')
        in  1 (mpf('1.000000000931322574615479'),) {}
        out 1 mpf('2.718281830990642467550102')
        mpf('2.71828')

    To disable either the input or the output handler, you may
    pass *None* as argument.

    Custom input and output handlers may be used e.g. to store
    results for later analysis::

        >>> mp.dps = 15
        >>> input = []
        >>> output = []
        >>> findroot(monitor(sin, input.append, output.append), 3.0)
        mpf('3.141592653589793')
        >>> len(input)  # Count number of evaluations
        9
        >>> print(input[3])
        ((mpf('3.1415076583334067'),), {})
        >>> print(output[3])
        8.499525628434076e-05
        >>> print(input[4])
        ((mpf('3.1415928201669123'),), {})
        >>> print(output[4])
        -1.6657711898533127e-07

    """
    if not input:
        input = lambda v: None
    elif input == 'print':
        incount = [0]
        def input(value):
            args, kwargs = value
            print("in  %s %r %r" % (incount[0], args, kwargs))
            incount[0] += 1
    if not output:
        output = lambda v: None
    elif output == 'print':
        outcount = [0]
        def output(value):
            print("out %s %r" % (outcount[0], value))
            outcount[0] += 1
    def f_monitored(*args, **kwargs):
        input((args, kwargs))
        v = f(*args, **kwargs)
        output(v)
        return v
    return f_monitored

def timing(f, *args, **kwargs):
    """
    Returns time elapsed for evaluating ``f()``. Optionally arguments
    may be passed to time the execution of ``f(*args, **kwargs)``.

    If the first call is very quick, ``f`` is called
    repeatedly and the best time is returned.
    """
    once = kwargs.get('once')
    if 'once' in kwargs:
        del kwargs['once']
    if args or kwargs:
        if len(args) == 1 and not kwargs:
            arg = args[0]
            g = lambda: f(arg)
        else:
            g = lambda: f(*args, **kwargs)
    else:
        g = f
    from timeit import default_timer as clock
    t1=clock(); v=g(); t2=clock(); t=t2-t1
    if t > 0.05 or once:
        return t
    for i in range(3):
        t1=clock()
        # Evaluate multiple times because the timer function
        # has a significant overhead
        g();g();g();g();g();g();g();g();g();g()
        t2=clock()
        t=min(t,(t2-t1)/10)
    return t
