import ast
import io
import tokenize


class IntegerDivisionWrapper(ast.NodeTransformer):
    """Wrap all int divisions in a call to :class:`~fractions.Fraction`."""

    def visit_BinOp(self, node):
        def is_integer(x):
            if isinstance(x, ast.Constant) and isinstance(x.value, int):
                return True
            if isinstance(x, ast.UnaryOp) and isinstance(x.op, (ast.USub,
                                                                ast.UAdd)):
                return is_integer(x.operand)
            if isinstance(x, ast.BinOp) and isinstance(x.op, (ast.Add,
                                                              ast.Sub,
                                                              ast.Mult,
                                                              ast.Pow)):
                return is_integer(x.left) and is_integer(x.right)
            return False

        if isinstance(node.op, ast.Div) and all(map(is_integer,
                                                    [node.left, node.right])):
            return ast.Call(ast.Name('Fraction', ast.Load()),
                            [node.left, node.right], [])
        return self.generic_visit(node)


def wrap_float_literals(lines):
    """Wraps all float/complex literals with mpmath classes."""
    result = []
    source = "\n".join(lines)
    g = tokenize.tokenize(io.BytesIO(source.encode()).readline)
    for toknum, tokval, *_ in g:
        if toknum == tokenize.NUMBER:
            if any(_ in tokval for _ in ['j', 'J']):
                tokval = f"mpc(imag=mpf('{tokval[:-1]}'))"
            elif any(_ in tokval for _ in ['.', 'e', 'E']):
                tokval = f"mpf('{tokval}')"
        result.append((toknum, tokval, *_))
    return tokenize.untokenize(result).decode().splitlines()
