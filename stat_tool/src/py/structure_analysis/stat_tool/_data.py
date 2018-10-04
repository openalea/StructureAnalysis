""" Low level wrapping of C++ functions into Python.

This module superseed data.py.
"""

from functools import wraps

from stat_tool import __stat_tool
import stat_tool.__stat_tool.stat_tool as cst


def wrapper(f):
    @wraps(f)
    def ascii_read(filename, *args, **kwds):
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        data = f(error, filename, *args, **kwds)
        if not data:
            raise Exception(str(error))
        return data
    return ascii_read


ascii_read_classes = [cst.DiscreteParametricModel]

for klass in ascii_read_classes:
    klass.ascii_read = staticmethod(wrapper(klass.ascii_read))

