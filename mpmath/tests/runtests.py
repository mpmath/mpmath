import py

# Makes test run over 50% faster
try:
    import psyco
    psyco.full()
except ImportError:
    pass

py.test.cmdline.main() 
