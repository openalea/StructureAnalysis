from functools import wraps

import __stat_tool
import _stat_tool

def wrapper(f):
    @wraps(f)
    def f_comparison(self, other):
        sstream = __stat_tool.std.Ostringstream(__stat_tool.std.ios_openmode.S__OUT)
        f(self, sstream, other)
        return sstream.str()
    return f_comparison

__stat_tool.stat_tool.FrequencyDistribution.f_comparison = wrapper(__stat_tool.stat_tool.FrequencyDistribution.f__comparison)
del __stat_tool.stat_tool.FrequencyDistribution.f__comparison

def wrapper(f):
    @wraps(f)
    def t_comparison(self, other):
        sstream = __stat_tool.std.Ostringstream(__stat_tool.std.ios_openmode.S__OUT)
        f(self, sstream, other)
        return sstream.str()
    return t_comparison

__stat_tool.stat_tool.FrequencyDistribution.t_comparison = wrapper(__stat_tool.stat_tool.FrequencyDistribution.t_comparison)

def wrapper(f):
    @wraps(f)
    def w_comparison(self, other):
        sstream = __stat_tool.std.Ostringstream(__stat_tool.std.ios_openmode.S__OUT)
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.NB_ERROR)
        f(self, error, sstream, other)
        if error.get_nb_error() > 0:
            raise Exception(str(error))
        return sstream.str()
    return w_comparison

__stat_tool.stat_tool.FrequencyDistribution.w_comparison = wrapper(__stat_tool.stat_tool.FrequencyDistribution.wilcoxon_mann_whitney_comparison)
del __stat_tool.stat_tool.FrequencyDistribution.wilcoxon_mann_whitney_comparison

def wrapper(f):
    @wraps(f)
    def discrete_mixture_estimation(self, imixt, min_inf_bound=0, mixt_flag=True, component_flag=True, weight_step=.1, **kwargs):
        pass
    return discrete_mixture_estimation


#__stat_tool.stat_tool.DiscreteDistributionData.ascii_read = staticmethod(__stat_tool.stat_tool.ascii_read_d00b270bd3835f568531951506bb4410)
#del __stat_tool.stat_tool.ascii_read_d00b270bd3835f568531951506bb4410

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