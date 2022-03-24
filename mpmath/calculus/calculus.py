class CalculusMethods(object):
    def __init__(self, *args, **kwargs):
        pass

def defun(f):
    setattr(CalculusMethods, f.__name__, f)
    return f
