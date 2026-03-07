"""Tests for demo scripts."""

import os
import subprocess
import sys
import time

import pexpect
import pytest


class Console(pexpect.spawn):
    """Spawned console for testing."""

    def __init__(self, command, timeout=60, _dumb=True):
        env = os.environ.copy()
        if _dumb:
            env['TERM'] = 'dumb'
        else:
            env['TERM'] = 'xterm'
            env['NO_COLOR'] = '1'
        super().__init__(command, timeout=timeout, encoding='utf-8', env=env)

    def __del__(self):
        self.send('exit()\r\n')
        time.sleep(10)  # a delay to allow coverage finish work
        if self.isalive():
            self.terminate(force=True)


# TODO: how to test plots? // mandelbrot.py and plotting.py


def test_manydigits():
    expected = r"""
This script prints answers to a selection of the "Many Digits"
competition problems: http://www.cs.ru.nl/~milad/manydigits/problems.php

The output for each problem is the first 100 digits after the
decimal point in the result.

C01: sin(tan(cos(1)))
56451092986195980582768640645029648577648661582588
56955552147245934844803576138875921296745208522197

C02: sqrt(e/pi)
93019136710263285866812462363333155602971092070428
87264450006489855422345460234483872155723942699765

C03: sin((e+1)^3)
90949524105726624718554721945217426889396524221380
80108799599078079083693175099387713504636663839042

C04: exp(pi*sqrt(2011))
08911292681099318912549002226654964403231616008375
14260187657441716605755144354088871641544234358651

C05: exp(exp(exp(1/2)))
33130360854569351505757451265398380886369247851475
92794392700131812592190818654155341658216570329325

C06: arctanh(1-arctanh(1-arctanh(1-arctanh(1/pi))))
12376761044118329658639748452701440281087636723733
55412845934779398491016984592299074199915669907895

C07: pi^1000
96790874439619754260235142488458363174182234378720
67532446047250097144332075967536835025898399733192

C08: sin(6^(6^6))
95395374345732063524921114340552534258118576365118
22065161716596988369691845451204872928519972839961

C09: sin(10*arctan(tanh(pi*(2011^(1/2))/3)))
99999999999999999999999999999999999999999999999999
99999999999999999999999999999868216408727535391618

C10: (7+2^(1/5)-5*(8^(1/5)))^(1/3) + 4^(1/5)-2^(1/5)
00000000000000000000000000000000000000000000000000
00000000000000000000000000000000000000000000000000

C11: tan(2^(1/2))+arctanh(sin(1))
56031033792570862486989423169964262718414115287379
65510969436882273871745968195963502918253580384966

C12: arcsin(1/e^2) + arcsinh(e^2)
83344680806041761874543293615785770019293386147122
63906848335142800750122119140978807925425237483497

C17: S= -4*Zeta(2) - 2*Zeta(3) + 4*Zeta(2)*Zeta(3) + 2*Zeta(5)
99922283776383000876193574924756988603699551613617
09442048984358627610229735501242221963535035597647

C18: Catalan G = Sum{i=0}{\infty}(-1)^i/(2i+1)^2
91596559417721901505460351493238411077414937428167
21342664981196217630197762547694793565129261151062

C21: Equation exp(cos(x)) = x
30296400121601255253211430697335802538621997810467
85962942111799929657676507417868401302803638230948

C22: J = integral(sin(sin(sin(x)))), x=0..1
40783902635001567262733691845249456720742376991339
01533400692321748591761662552762179981626145798049

"""
    result = subprocess.run([f'{sys.executable}',
                             'demo/manydigits.py'],
                            capture_output=True, text=True)
    assert result.stdout == expected


def test_pidigits():
    c = Console(f'{sys.executable} demo/pidigits.py')
    assert c.expect_exact('> ') == 0
    assert c.send('10\n') == 3
    assert c.expect_exact('> ') == 0
    assert c.send('100\n') == 4
    assert c.expect_exact('> ') == 0
    assert c.send('\n') == 1
    assert c.expect('5820974944 5923078164 0628620899 '
                    '8628034825 3421170679 : 100') == 0


def test_sofa():
    result = subprocess.run([f'{sys.executable}',
                             'demo/sofa.py'],
                            capture_output=True, text=True)
    assert result.stdout == '2.2195316688719674255462841007968\n'


def test_taylor():
    c = Console(f'{sys.executable} demo/taylor.py')
    assert c.expect_exact('Enter the value of x (e.g. 3.5): ') == 0
    assert c.send('1\n') == 2
    assert c.expect_exact('Enter the number of terms n (e.g. 10): ') == 0
    assert c.send('10\n') == 3
    assert c.expect_exact('[2.7182818011463827368, 2.7182818011463862895]') == 0
