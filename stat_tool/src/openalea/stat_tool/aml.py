def load_ipython_extension(ipython):

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
        return cls(args)

    stat_tool.Merge = Merge
    del Merge

    def Plot(obj, *args, **kwargs):
        return obj.plot(*args, **kwargs)


def unload_ipython_extension(ipython):

    from openalea.stat_tool import stat_tool

    del stat_tool.Histogram

