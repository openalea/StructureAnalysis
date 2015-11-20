import sys

DISABLE_PLOT = False


# !!! Do not plot in nosetests !!!
# buildbot cannot close the windwos popped up by the method/function Plot/plot
# So, we test if the command "python setup.py nosetests" has been used.
# Still, using nosetests executable, windows should pop up.
if("nosetests" in sys.argv):
    DISABLE_PLOT = True

class plotter(object):
    """ Abstract base class for all plotter """
    def _init__(self):
        pass

    def plot(self, obj, title, groups=None, *args, **kargs):
        """ Plot obj with title """

        raise NotImplementedError()


class fakeplot(plotter):

    def __init__(self):
        plotter.__init__(self)

    def plot(self, obj, title, groups=None, *args, **kargs):
        """ Plot obj with title """
        return


class gnuplot(plotter):
    """ GNUPlot implementation """

    def __init__(self):
        """ Initialize GnuPlot """
        plotter.__init__(self)
        import Gnuplot
        self.session = Gnuplot.Gnuplot()

    def plot(self, plotable, title, groups=None, *args, **kargs):
        """
         Plot a plotable with title
         groups : list of group (int) to plot
        """
        import Gnuplot

        multiset = plotable
        g = self.session

        # Title & border
        #multiset.border
        # nb subplot
        _nbx = len(multiset)

        # For each subplot
        for _i, multiplot in enumerate(multiset):
            # Group filter
            if(groups and multiplot.group not in groups):
                continue
            g.title(multiplot.title)
            #yoffset = i * plotsize

            #g('set origin 0.0, %f'%(yoffset))
            #g('set size 1.0, %f'%(plotsize))

            # Labels
            g.xlabel(multiplot.xlabel)
            g.ylabel(multiplot.ylabel)

            # List of argument for the plot function
            plot_list = []
            for singleplot in multiplot:

                style = singleplot.style
                legend = singleplot.legend
                _color = singleplot.color

                x = []
                y = []
                for pt in singleplot:
                    x.append(pt.first)
                    y.append(pt.second)

                p = Gnuplot.Data(x, y)
                if style:
                    #todo: check that this option works.
                    p.set_option(with_=style)

                if legend:
                    p.set_option(title=legend)
                plot_list.append(p)

            # Range
            _xrange = multiplot.xrange
            _yrange = multiplot.yrange
            if(_xrange.min != _xrange.max):
                g('set xrange[%f:%f]'%(_xrange.min, _xrange.max))
            if(_yrange.min != _yrange.max):
                g('set yrange[%f:%f]'%(_yrange.min, _yrange.max))

            # Tics
            if(multiplot.xtics > 0):
                g('set xtics 0, %f'%(multiplot.xtics))
            if(multiplot.ytics > 0):
                g('set ytics 0, %f'%(multiplot.ytics))


            g.plot(*plot_list)
            raw_input("Press Enter to continue")


PLOTTER = None


def set_plotter(plot):
    global PLOTTER
    PLOTTER = plot


def get_plotter():
    """
    Plotter factory
    Return a plotter object (matplotlib or gnuplot)
    If none is available, raise an ImportError exception
    """
    global DISABLE_PLOT
    if DISABLE_PLOT:
        return fakeplot()

    # Try user define PLOTTER
    global PLOTTER
    if PLOTTER:
        return PLOTTER

    # Try to import a plotter
    else:
        try:
            _plotter = mplotlib()

        except ImportError:
            _plotter = gnuplot()

        set_plotter(_plotter)
        return _plotter
