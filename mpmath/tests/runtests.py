#!python

"""
python runtests.py -py
  Use py.test to run tests (more useful for debugging)

python runtests.py -psyco
  Enable psyco to make tests run about 50% faster

python runtests.py -profile
  Generate profile stats (this is much slower)

"""

import sys

if "-psyco" in sys.argv:
    sys.argv.remove('-psyco')
    import psyco
    psyco.full()

profile = False
if "-profile" in sys.argv:
    sys.argv.remove('-profile')
    profile = True

def testit():
    if "-py" in sys.argv:
        sys.argv.remove('-py')
        import py
        py.test.cmdline.main()
    else:
        import glob
        import os.path
        from time import clock
        modules = []
        for f in glob.glob("test*.py"):
            name = os.path.splitext(os.path.basename(f))[0]
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

if profile:
    import cProfile
    cProfile.run("testit()", sort=2)
else:
    testit()
