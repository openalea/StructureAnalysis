from functools import wraps

from openalea.stat_tool import __stat_tool

__stat_tool.stat_tool.DiscreteDistributionData.ascii_read = staticmethod(__stat_tool.stat_tool.ascii_read_d00b270bd3835f568531951506bb4410)
del __stat_tool.stat_tool.ascii_read_d00b270bd3835f568531951506bb4410

def wrapper(f):
    @wraps(f)
    def ascii_read(filename):
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.NB_ERROR)
        data = f(error, filename)
        if not data:
            raise Exception(str(error))
        return data
    return ascii_read

__stat_tool.stat_tool.DiscreteDistributionData.ascii_read = staticmethod(wrapper(__stat_tool.stat_tool.DiscreteDistributionData.ascii_read))
del wrapper
