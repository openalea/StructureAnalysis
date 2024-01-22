def load_ipython_extension(ipython):

    from multipledispatch import dispatch

    import stat_tool
    from stat_tool import __stat_tool

    def Histogram(filename):
        """
        """
        return __stat_tool.stat_tool.DiscreteDistributionData.ascii_read(filename)

    stat_tool.Histogram = Histogram
    del Histogram

    def Merge(*args):
        if not all(hasattr(arg, 'merge') for arg in args):
            raise AttributeError('no `merge` method found')
        return args[0].merge(*args[1:])

    stat_tool.Merge = Merge
    del Merge

    def Plot(obj, *args, **kwargs):
        return obj.plot(*args, **kwargs)

    @dispatch(str, __stat_tool.stat_tool.DiscreteDistributionData, __stat_tool.stat_tool.DiscreteDistributionData)
    def ComparisonTest(test, lhs, rhs):
        if test == 'F':
            test = lhs.f_comparison(rhs)
        elif test == 'T':
            test = lhs.t_comparison(rhs)
        elif test == 'W':
            test = lhs.w_comparison(rhs)
        else:
            raise ValueError('\'test\' parameter')
        print(test)

    stat_tool.ComparisonTest = ComparisonTest

    def MixtureEstimation(data, *args, **kwargs):
        if isinstance(data, __stat_tool.stat_tool.DiscreteDistributionData):
            return data.discrete_mixture_estimation(*args, **kwargs)
        else:
            raise NotImplementedError('For a `' + data.__class__.__name__ + '` instance')

    stat_tool.MixtureEstimation = MixtureEstimation

def unload_ipython_extension(ipython):

    import stat_tool

    del stat_tool.Histogram
