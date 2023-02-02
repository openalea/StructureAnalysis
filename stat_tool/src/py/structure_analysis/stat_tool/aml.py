def load_ipython_extension(ipython):

    from multipledispatch import dispatch

    from openalea.stat_tool import stat_tool
    from openalea.stat_tool import __stat_tool

    def Histogram(filename):
        """
        """
        return __stat_tool.stat_tool.DiscreteDistributionData.ascii_read(filename)

    stat_tool.Histogram = Histogram
    del Histogram

    def Merge(*args):
        cls = args[0].__class__
        if not all(isinstance(arg, cls) for arg in args):
            raise TypeError('Cannot merge objects of different types')
        return cls(args)

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
        print test

    stat_tool.ComparisonTest = ComparisonTest


def unload_ipython_extension(ipython):

    from openalea.stat_tool import stat_tool

    del stat_tool.Histogram