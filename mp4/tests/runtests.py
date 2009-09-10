#!/usr/bin/env python

"""
python runtests.py -py
  Use py.test to run tests (more useful for debugging)

python runtests.py -psyco
  Enable psyco to make tests run about 50% faster

python runtests.py -coverage
  Generate test coverage report. Statistics are written to /tmp

python runtests.py -profile
  Generate profile stats (this is much slower)

python runtests.py -nogmpy
  Run tests without using GMPY even if it exists

python runtests.py -strict
  Enforce extra tests in normalize()

python runtests.py -local
  Insert '../..' at the beginning of sys.path to use local mpmath

Additional arguments are used to filter the tests to run. Only files that have
one of the arguments in their name are executed.

"""

import sys, os

if "-psyco" in sys.argv:
    sys.argv.remove('-psyco')
    import psyco
    psyco.full()

profile = False
if "-profile" in sys.argv:
    sys.argv.remove('-profile')
    profile = True

coverage = False
if "-coverage" in sys.argv:
    sys.argv.remove('-coverage')
    coverage = True

if "-nogmpy" in sys.argv:
    sys.argv.remove('-nogmpy')
    os.environ['MPMATH_NOGMPY'] = 'Y'

if "-strict" in sys.argv:
    sys.argv.remove('-strict')
    os.environ['MPMATH_STRICT'] = 'Y'

if "-local" in sys.argv:
    sys.argv.remove('-local')
    importdir = os.path.abspath(os.path.join(os.path.dirname(sys.argv[0]),
                                             '../..'))
else:
    importdir = ''

# TODO: add a flag for this
testdir = ''

def testit(importdir='', testdir=''):
    """Run all tests in testdir while importing from importdir."""
    if importdir:
        sys.path.insert(1, importdir)
    if testdir:
        sys.path.insert(1, testdir)
    if "-py" in sys.argv:
        sys.argv.remove('-py')
        import py
        py.test.cmdline.main()
    else:
        import glob
        import os.path
        from timeit import default_timer as clock
        modules = []
        args = sys.argv[1:]
        # search for tests in directory of this file if not otherwise specified
        if not testdir:
            pattern = os.path.dirname(sys.argv[0])
        else:
            pattern = testdir
        if pattern:
            pattern += '/'
        pattern += 'test*.py'
        # look for tests (respecting specified filter)
        for f in glob.glob(pattern):
            name = os.path.splitext(os.path.basename(f))[0]
            # If run as a script, only run tests given as args, if any are given
            if args and __name__ == "__main__":
                ok = False
                for arg in args:
                    if arg in name:
                        ok = True
                        break
                if not ok:
                    continue
            module = __import__(name)
            priority = module.__dict__.get('priority', 100)
            if priority == 666:
                modules = [[priority, name, module]]
                break
            modules.append([priority, name, module])
        # execute tests
        modules.sort()
        tstart = clock()
        for priority, name, module in modules:
            print name
            for f in sorted(module.__dict__.keys()):
                if f.startswith('test_'):
                    if coverage and ('numpy' in f):
                        continue
                    print "   ", f[5:].ljust(25),
                    t1 = clock()
                    module.__dict__[f]()
                    t2 = clock()
                    print "ok", "      ", ("%.7f" % (t2-t1)), "s"
        tend = clock()
        print
        print "finished tests in", ("%.2f" % (tend-tstart)), "seconds"
        # clean sys.path
        if importdir:
            sys.path.remove(importdir)
        if testdir:
            sys.path.remove(testdir)

if __name__ == '__main__':
    if profile:
        import cProfile
        cProfile.run("testit('%s', '%s')" % (importdir, testdir), sort=1)
    elif coverage:
        import trace
        tracer = trace.Trace(ignoredirs=[sys.prefix, sys.exec_prefix],
            trace=0, count=1)
        tracer.run('testit(importdir, testdir)')
        r = tracer.results()
        r.write_results(show_missing=True, summary=True, coverdir="/tmp")
    else:
        testit(importdir, testdir)

