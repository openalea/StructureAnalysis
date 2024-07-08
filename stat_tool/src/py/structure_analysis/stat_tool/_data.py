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


_classes = [cst.DiscreteParametricModel, cst.DiscreteMixture]

for klass in _classes:
    klass.ascii_read = staticmethod(wrapper(klass.ascii_read))

del wrapper

def wrapper(f):
    @wraps(f)
    def simulation(self, size, *args, **kwds):
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        data = f(self, error, size, *args, **kwds)
        if not data:
            raise Exception(str(error))
        return data
    return simulation

_classes = [cst.DiscreteParametricModel, cst.DiscreteMixture]
for klass in _classes:
    klass.simulation = wrapper(klass.simulation)

##############################################################################

def wrapper(f):
    @wraps(f)
    def parametric_estimation(self, *args, **kwds):
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        data = f(self, error,*args, **kwds)
        if not data:
            raise Exception(str(error))
        return data
    return parametric_estimation

_classes = [cst.FrequencyDistribution]
for klass in _classes:
    klass.parametric_estimation = wrapper(klass.parametric_estimation)



def wrapper(f):
    @wraps(f)
    def __str__(self):
        data = self.ascii_write(False)
        return data
    return __str__

def __ascii_w(module):
    for obj in module.__dict__.itervalues():
        if hasattr(obj, 'ascii_write'):
            obj.__str__ = wrapper(getattr(obj, 'ascii_write'))
_classes =__ascii_w(cst)
