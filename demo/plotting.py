"""
Function plotting demo.
"""

from mpmath import *

def main():
    print """
    Simple function plotting. You can enter one or several
    formulas, in ordinary Python syntax and using the mpmath
    function library. The variable is 'x'. So for example
    the input "sin(x/2)" (without quotation marks) defines
    a valid function.
    """
    functions = []
    for i in xrange(10):
        if i == 0:
            s = raw_input('Enter a function: ')
        else:
            s = raw_input('Enter another function (optional): ')
        if not s:
            print
            break
        f = eval("lambda x: " + s)
        functions.append(f)
        print "Added f(x) = " + s
        print
    xlim = raw_input('Enter xmin, xmax (optional): ')
    if xlim:
        xlim = eval(xlim)
    else:
        xlim = [-5, 5]
    print "Plotting..."
    plot(functions, xlim=xlim)

main()