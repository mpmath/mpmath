"""
Python shell for mpmath.

This is just a normal Python shell (IPython shell if you have the
IPython package installed), that adds default imports and run
some initialization code.
"""

import argparse
import ast
import atexit
import code
import os
import readline
import rlcompleter
import sys

from mpmath import __version__
from mpmath._interactive import (IntegerDivisionWrapper,
                                 wrap_float_literals)


__all__ = ()


parser = argparse.ArgumentParser(description=__doc__,
                                 prog='python -m mpmath')
parser.add_argument('--no-wrap-division',
                    help="Don't wrap integer divisions with Fraction",
                    action='store_true')
parser.add_argument('--no-ipython', help="Don't use IPython",
                    action='store_true')
parser.add_argument('--no-wrap-floats',
                    help="Don't wrap float/complex literals",
                    action='store_true')
parser.add_argument('-V', '--version',
                    help='Print the mpmath version and exit',
                    action='store_true')
parser.add_argument('--prec', type=int,
                    help='Set default mpmath precision')
parser.add_argument('--pretty', help='Enable pretty-printing',
                    action='store_true')


def main():
    args, ipython_args = parser.parse_known_args()

    if args.version:
        print(__version__)
        sys.exit(0)

    lines = ['from mpmath import *',
             'from fractions import Fraction']

    if args.prec:
        lines.append(f'mp.prec = {args.prec}')
    if args.pretty:
        lines.append('mp.pretty = True')

    try:
        import IPython
        import traitlets
    except ImportError:
        args.no_ipython = True

    if not args.no_ipython:
        config = traitlets.config.loader.Config()
        shell = config.InteractiveShell
        ast_transformers = shell.ast_transformers
        if not args.no_wrap_division:
            ast_transformers.append(IntegerDivisionWrapper())
        shell.confirm_exit = False
        config.TerminalIPythonApp.display_banner = False
        config.TerminalInteractiveShell.autoformatter = None

        app = IPython.terminal.ipapp.TerminalIPythonApp.instance(config=config)
        app.initialize(ipython_args)
        shell = app.shell
        for l in lines:
            shell.run_cell(l, silent=True)
        if not args.no_wrap_floats:
            shell.run_cell('from mpmath._interactive import wrap_float_literals')
            shell.run_cell('ip = get_ipython()')
            shell.run_cell('ip.input_transformers_post.append(wrap_float_literals)')
            shell.run_cell('del ip')
        app.start()
    else:
        ast_transformers = []
        source_transformers = []
        ns = {}

        if not args.no_wrap_division:
            ast_transformers.append(IntegerDivisionWrapper())
        if not args.no_wrap_floats:
            source_transformers.append(wrap_float_literals)

        class MpmathConsole(code.InteractiveConsole):
            """An interactive console with readline support."""

            def __init__(self, ast_transformers=[],
                         source_transformers=[], **kwargs):
                super().__init__(**kwargs)

                readline.set_completer(rlcompleter.Completer(ns).complete)
                readline.parse_and_bind('tab: complete')

                history = os.path.expanduser('~/.python_history')
                readline.read_history_file(history)
                atexit.register(readline.write_history_file, history)
                self.ast_transformers = ast_transformers
                self.source_transformers = source_transformers

            def runsource(self, source, filename='<input>', symbol='single'):
                if not source:
                    return True

                for t in self.source_transformers:
                    source = '\n'.join(t(source.splitlines()))

                if self.ast_transformers:
                    tree = ast.parse(source, mode=symbol)

                    for t in self.ast_transformers:
                        tree = t.visit(tree)
                    ast.fix_missing_locations(tree)

                    code_obj = compile(tree, filename, mode=symbol)
                    try:
                        self.runcode(code_obj)
                    except SystemExit:
                        os.exit(0)
                    return False

                return super().runsource(source, filename=filename, symbol=symbol)

        c = MpmathConsole(ast_transformers=ast_transformers,
                          source_transformers=source_transformers, locals=ns)

        for l in lines:
            c.push(l)
        c.interact('', '')


if __name__ == '__main__':
    main()
