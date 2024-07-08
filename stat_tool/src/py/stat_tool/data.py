from functools import wraps

from . import _stat_tool
from . import __stat_tool

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
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        f(self, error, sstream, other)
        if error.get_nb_error() > 0:
            raise Exception(str(error))
        return sstream.str()
    return w_comparison

__stat_tool.stat_tool.FrequencyDistribution.w_comparison = wrapper(__stat_tool.stat_tool.FrequencyDistribution.wilcoxon_mann_whitney_comparison)
del __stat_tool.stat_tool.FrequencyDistribution.wilcoxon_mann_whitney_comparison

def wrapper(f):
    @wraps(f)
    def merge(self, *args):
        for arg in args:
            if not isinstance(arg, __stat_tool.stat_tool.FrequencyDistribution):
                raise TypeError('cannot merge `' + arg.__class__.__name__ + '` instance to `' + self.__class__.__name__ + '` instance')
        return f(self, len(args), args)
    return merge

__stat_tool.stat_tool.FrequencyDistribution.merge = wrapper(__stat_tool.stat_tool.FrequencyDistribution.merge)
del wrapper

def wrapper(f):
    @wraps(f)
    def discrete_mixture_estimation(self, max_nb_component, ident, min_inf_bound=0,
                                    mixt_flag=True, component_flag=True,
                                    model_selection_criterion="BIC", weight_step=.1, **kwargs):
        if isinstance(ident, basestring):
            ident = [ident] * max_nb_component
        if isinstance(ident, (list, tuple)):
            ident = [__stat_tool.stat_tool.discrete_parametric.names[ident] for ident in ident]
        model_selection_criterion = __stat_tool.stat_tool.model_selection_criterion.names[model_selection_criterion]
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        mixt = f(self,
                 error,
                 kwargs.pop('display', False),
                 kwargs.pop('min_nb_component', 1),
                 max_nb_component,
                 ident,
                 min_inf_bound,
                 mixt_flag,
                 component_flag,
                 model_selection_criterion,
                 weight_step)
        if error.get_nb_error() > 0:
            raise Exception(str(error))
        return mixt
    return discrete_mixture_estimation

__stat_tool.stat_tool.FrequencyDistribution.discrete_mixture_estimation = wrapper(__stat_tool.stat_tool.FrequencyDistribution.discrete_mixture_estimation)

def wrapper(f):
    @wraps(f)
    def ascii_read(filename):
        error = __stat_tool.stat_tool.StatError(__stat_tool.stat_tool.nb_error)
        data = f(error, filename)
        if not data:
            raise Exception(str(error))
        return data
    return ascii_read

__stat_tool.stat_tool.DiscreteDistributionData.ascii_read = staticmethod(wrapper(__stat_tool.stat_tool.DiscreteDistributionData.ascii_read))
del wrapper
