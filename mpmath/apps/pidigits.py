"""
Calculate digits of pi using AGM iteration. The implementation relies on
mpmath.lib's fixed-point square root code to do most of the work.

This module can be run interactively with

    python pidigits.py

or from the interpreter using

    >>> from mpmath.apps.pidigits import *
    >>> interactive()

"""

import sys
import math
from time import clock

from mpmath.lib import bin_to_radix, numeral, sqrt_fixed2

def agm_status(prec, step, adiff, verbose_base):
    logdiff = math.log(max(1, adiff), verbose_base)
    digits = int(prec/math.log(verbose_base,2) - logdiff)
    print "  iteration", step, ("(accuracy ~= %i base-%i digits)" % \
       (digits, verbose_base))

def pi_agm(prec, verbose=False, verbose_base=10):
    """
    Compute floor(pi * 2**prec) as a big integer using the Brent-
    Salamin algorithm based on the arithmetic-geometric mean.

    See for example Wikipedia (http://en.wikipedia.org/wiki/Brent-
    Salamin_algorithm) or "Pi and the AGM" by Jonathan and Peter
    Borwein (Wiley, 1987). The algorithm (as stated in the Wikipedia
    article) consists of setting

      a_0 = 1
      b_0 = 1/sqrt(2)
      t_0 = 1/4
      p_0 = 1

    and computing

      a_{n+1} = (a_n + b_n)/2
      b_{n+1} = sqrt(a_n * b_n)
      t_{n+1} = t_n - p_n*(a_n - a_{n+1})**2
      p_{n+1} = 2*p_n

    for n = 0, 1, 2, 3, ..., after which the approximation is given by
    pi ~= (a_n + b_n)**2 / (4*t_n). Each step roughly doubles the
    number of correct digits.
    """
    extraprec = 50
    prec += extraprec

    # Initialial values. a, b and t are fixed-point numbers
    a = 1 << prec
    b = sqrt_fixed2(a >> 1, prec)
    t = a >> 2
    p = 1

    step = 1
    while 1:
        an = (a + b) >> 1
        adiff = a - an
        if verbose:
            agm_status(prec, step, adiff, verbose_base)
        # No change in a
        if p > 16 and abs(adiff) < 1000:
            break
        prod = (a * b) >> prec
        b = sqrt_fixed2(prod, prec)
        t = t - p*((adiff**2) >> prec)
        p = 2*p
        a = an
        step += 1
    if verbose:
        print "  final division"
    pi = ((((a+b)**2) >> 2) // t)
    return pi >> extraprec

def display_fraction(digits, skip=0, colwidth=10, columns=5):
    perline = colwidth * columns
    printed = 0
    for linecount in range((len(digits)-skip) // (colwidth * columns)):
        line = digits[skip+linecount*perline:skip+(linecount+1)*perline]
        for i in range(columns):
            print line[i*colwidth : (i+1)*colwidth],
        print ":", (linecount+1)*perline
        if (linecount+1) % 10 == 0:
            print
        printed += colwidth*columns
    rem = (len(digits)-skip) % (colwidth * columns)
    if rem:
        buf = digits[-rem:]
        s = ""
        for i in range(columns):
            s += buf[:colwidth].ljust(colwidth+1, " ")
            buf = buf[colwidth:]
        print s + ":", printed + colwidth*columns

def calculateit(base, n, tofile):
    intpart = numeral(3, base)
    skip = 1

    prec = int(n*math.log(base,2))+10

    print "Step 1 of 2: calculating binary value..."
    t = clock()
    a = pi_agm(prec, verbose=True, verbose_base=base)
    step1_time = clock() - t

    print "Step 2 of 2: converting to specified base..."
    t = clock()
    d = bin_to_radix(a, prec, base, n)
    d = numeral(d, base, n)
    step2_time = clock() - t

    print "\nWriting output...\n"

    if tofile:
        out_ = sys.stdout
        sys.stdout = tofile
    print "%i base-%i digits of pi:\n" % (n, base)
    print intpart, ".\n"

    display_fraction(d, skip, colwidth=10, columns=5)
    if tofile:
        sys.stdout = out_
    print "\nFinished in %f seconds (%f calc, %f convert)" % \
        ((step1_time + step2_time), step1_time, step2_time)

def interactive():
    print "Compute digits of pi with mpmath\n"
    base = input("Which base? (2-36, 10 for decimal) \n> ")
    digits = input("How many digits? (enter a big number, say, 10000)\n> ")
    tofile = raw_input("Output to file? (enter a filename, or just press " \
        "enter\nto print directly to the screen) \n> ")
    if tofile:
        tofile = open(tofile, "w")

    calculateit(base, digits, tofile)
    raw_input("\nPress enter to close this script.")

if __name__ == "__main__":
    interactive()
