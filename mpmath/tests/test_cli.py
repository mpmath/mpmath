"""Tests for the Command-Line Interface."""

import platform
import sys
import time

import pexpect
import pytest


if platform.python_implementation() == 'PyPy':
    pytest.skip("Don't run CLI tests on PyPy.",
                allow_module_level=True)


class Console(pexpect.spawn):
    """Spawned console for testing."""

    def __init__(self, command, timeout=60):
        super().__init__(command, timeout=timeout, encoding='utf-8')

    def __del__(self):
        self.send('exit()\r\n')
        time.sleep(10)  # a delay to allow coverage finish work
        if self.isalive():
            self.terminate(force=True)


def test_bare_console_no_bare_division():
    c = Console(f'{sys.executable} -m mpmath --no-ipython --no-wrap-floats')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1 + 2\r\n') == 7
    assert c.expect_exact('3\r\n>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('Fraction(1, 2)\r\n>>> ') == 0
    assert c.send('-1/2\r\n') == 6
    assert c.expect_exact('Fraction(-1, 2)\r\n>>> ') == 0
    assert c.send('2**3/7\r\n') == 8
    assert c.expect_exact('Fraction(8, 7)\r\n>>> ') == 0
    assert c.send('(3 + 5)/7\r\n') == 11
    assert c.expect_exact('Fraction(8, 7)\r\n>>> ') == 0
    assert c.send('(0.5 + 1)/2\r\n') == 13
    assert c.expect_exact('0.75\r\n>>> ') == 0


def test_bare_console_bare_division():
    c = Console(f'{sys.executable} -m mpmath --no-ipython --no-wrap-division '
                '--no-wrap-floats')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('0.5\r\n>>> ') == 0


def test_bare_console_without_ipython():
    try:
        import IPython
        del IPython
        pytest.skip('IPython is available')
    except ImportError:
        pass

    c = Console(f'{sys.executable} -m mpmath')

    assert c.expect_exact('>>> ') == 0
    assert c.send('1 + 2\r\n') == 7
    assert c.expect_exact('3\r\n>>> ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nFraction(1, 2)\r\n>>> ') == 0


def test_ipython_console_bare_division_noauto():
    pytest.importorskip('IPython')

    c = Console(f'{sys.executable} -m mpmath --simple-prompt --no-wrap-floats '
                "--no-wrap-division --colors 'NoColor' ")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('1/2\r\n') == 5
    assert c.expect_exact('\r\nOut[1]: 0.5\r\n\r\nIn [2]: ') == 0


def test_ipython_console_wrap_floats():
    pytest.importorskip('IPython')

    c = Console(f'{sys.executable} -m mpmath --simple-prompt --prec 100 '
                "--colors 'NoColor'")

    assert c.expect_exact('\r\nIn [1]: ') == 0
    assert c.send('10.9\r\n') == 6
    assert c.expect_exact("\r\nOut[1]: mpf('10.899999999999999999999999999995')\r\n\r\nIn [2]: ") == 0


def test_bare_console_wrap_floats():
    c = Console(f'{sys.executable} -m mpmath --simple-prompt --no-ipython --prec 100 '
                "--colors 'NoColor'")

    assert c.expect_exact('>>> ') == 0
    assert c.send("10.9\r\n") == 6
    assert c.expect_exact("mpf('10.899999999999999999999999999995')\r\n>>> ") == 0
    assert c.send("1+10.9j\r\n") == 9
    assert c.expect_exact("mpc(real='1.0', imag='10.899999999999999999999999999995')\r\n>>> ") == 0
    assert c.send('mpf(10.9)\r\n') == 11
    assert c.expect_exact("mpf('10.899999999999999999999999999995')\r\n>>> ") == 0


def test_bare_console_pretty():
    c = Console(f'{sys.executable} -m mpmath --simple-prompt --no-ipython --prec 100 '
                "--colors 'NoColor' --pretty")

    assert c.expect_exact('>>> ') == 0
    assert c.send("10.9\r\n") == 6
    assert c.expect_exact("10.9\r\n>>> ") == 0


def test_mpmath_version():
    c = Console(f'{sys.executable} -m mpmath --version')

    assert c.expect(pexpect.EOF) == 0
    assert c.before.startswith('1.')
