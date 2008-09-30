#!/usr/bin/env python

"""
python runtests.py -py
  Use py.test to run tests (more useful for debugging)

python runtests.py -psyco
  Enable psyco to make tests run about 50% faster

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

if "-nogmpy" in sys.argv:
    sys.argv.remove('-nogmpy')
    os.environ['MPMATH_NOGMPY'] = 'Y'

if "-strict" in sys.argv:
    sys.argv.remove('-strict')
    os.environ['MPMATH_STRICT'] = 'Y'

if "-local" in sys.argv:
    sys.argv.remove('-local')
    directory = '../..'
else:
    directory = ''

# TODO: make it possible to run it from another directory
def testit(directory=''):
    """Run all tests while importing from directory."""
    if directory:
        sys.path.insert(1, directory)
    if "-py" in sys.argv:
        sys.argv.remove('-py')
        import py
        py.test.cmdline.main()
    else:
        import glob
        import os.path
        from time import clock
        modules = []
        args = sys.argv[1:]
        for f in glob.glob("test*.py"):
            name = os.path.splitext(os.path.basename(f))[0]
            print name
            if args:
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
        modules.sort()
        tstart = clock()
        for priority, name, module in modules:
            print name
            for f in sorted(module.__dict__.keys()):
                if f.startswith('test_'):
                    print "   ", f[5:].ljust(25),
                    t1 = clock()
                    module.__dict__[f]()
                    t2 = clock()
                    print "ok", "      ", ("%.7f" % (t2-t1)), "s"
        tend = clock()
        print
        print "finished tests in", ("%.2f" % (tend-tstart)), "seconds"

if __name__ == '__main__':
    if profile:
        import cProfile
        cProfile.run("testit(%s)" % directory, sort=2)
    else:
        testit(directory)

