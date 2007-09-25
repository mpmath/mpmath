import sys

import math
from time import clock

from mpmath import *
from mpmath.lib import _pi_agm

def fpi(prec=STANDARD_PREC, rounding=ROUND_HALF_EVEN):
    """Compute a floating-point approximation of pi"""
    return normalize(pi_fixed(prec+5), -prec-5, prec, rounding)


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

def calculateit(approx, func, base, n, tofile):

    intpart = small_numeral(int(approx), base)
    if intpart == 0:
        skip = 0
    else:
        skip = len(intpart)

    prec = int(n*math.log(base,2))+10
    mpf.prec = prec

    print "Step 1 of 2: calculating binary value..."
    t = clock()
    a = func(prec=prec, verbose=True, verbose_base=base)
    a = mpf(normalize(a, -prec, prec))
    step1_time = clock() - t

    print "Step 2 of 2: converting to specified base..."
    t = clock()
    d = bin_to_radix(a.man, -a.exp, base, n)
    d = fixed_to_str(d, base, n)
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
    tofile = raw_input("Output to file? (enter a filename, or just press enter\nto print directly to the screen) \n> ")
    if tofile:
        tofile = open(tofile, "w")

    calculateit(3.14, _pi_agm, base, digits, tofile)
    raw_input("\nPress enter to close this script.")

interactive()
