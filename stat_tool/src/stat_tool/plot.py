#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Plot functions

.. topic:: plot.py summary

    A module that provides plotting functions for Gnuplot
    of Matplorlib.

    :Code status: mature
    :Documentation status: to be completed
    :Authors: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>, 
        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>

    :Revision: $Id: plot.py 13231 2012-11-27 15:48:09Z jbdurand $

.. inheritance-diagram:: openalea.stat_tool.plot

"""
__version__ = "$Id: plot.py 13231 2012-11-27 15:48:09Z jbdurand $"

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
                    x.append(pt.x)
                    y.append(pt.y)

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


class mplotlib(plotter):
    """matplotlib implementation of AML Plot


    This class defines an interface to matplolib plotting functions
    and options.



    """

    linestyles = ('-', '--', ':', '.')
    pointstyles = ('o', '^', 'x', '+', 's', 'v', '>', '<')
    colors = ('g', 'r', 'y', 'b', 'k', 'm', 'c')

    def __init__(self):
        """ Initialize matplotlib """
        plotter.__init__(self)
        import matplotlib
        #matplotlib.use('Qt4Agg')
        import pylab
        self.pylab = pylab
        self.matplotlib = matplotlib

    def plot(self, plotable, title, groups=None, *args, **kargs):
        """Plot a plotable with title

        :param plotable: a plotable instance from standard stat_tool objects such as :func:`Distribution`
        :param groups: list of group (int) to plot
        :param Show:
        :param nbcol:
        :param legend_size: 10
        :param legend_nbcol: 2
        :param legend_loc: best
        :param legend: True/False

        show=True by default will pop up the figure
        y_maxrange_ratio=1 multiply max y range by this value
        """

        show = kargs.get("Show", True)
        nbcol = kargs.get("nbcol", 2)
        fontsize = kargs.get("fontsize", 10)
        options_y_maxrange_ratio = kargs.get("y_maxrange_ratio", 1.1)
        line2d = {}
        line2d['linewidth'] = kargs.get("linewidth", 1)
        legend_size = kargs.get("legend_size", 10)
        legend_nbcol = kargs.get("legend_nbcol", 1)
        legend_loc = kargs.get("legend_loc", "best")
        legend_on = kargs.get("legend", True)
        #legend_size = kargs.get("legend_kwds", {})
        #Plot(seq1, ViewPoint="Data", nbcol=2, legend_kwds={'prop':{'size':9}})
        fig_id = kargs.get("FigureId", 1)

        pylab = self.pylab
        matplotlib = self.matplotlib
        multiset = plotable


        if len(multiset)==1:
            nbcol = 1
        # Title & border
        #if title: pylab.suptitle(title)
        #multiset.border

        # Configure figure
        f1 = pylab.figure(fig_id, figsize=(10, 10))
        f1.clf()

        f1.set_facecolor("w")

        # Count group
        group_size = {} # Group size
        group_index = {} # group counter
        for multiplot in multiset:
            g = multiplot.group
            try:
                group_size[g] += 1

            except KeyError:
                # init group
                group_size[g] = 1
                group_index[g] = 0

                f = pylab.figure(fig_id + g )
                f.set_facecolor("w")

        # nb subplot
        _nbx = len(multiset)

        if title:
            pylab.suptitle(title, fontsize=fontsize)
        # For each subplot
        for i, multiplot in enumerate(multiset):
            g = multiplot.group
            # Group filter
            if(groups and g not in groups):
                continue

            # Select window
            pylab.figure(fig_id + g )

            # Select Subplot
            pylab.subplot(group_size[g]%nbcol+group_size[g]/nbcol,
                          nbcol , group_index[g] + 1)
            group_index[g] += 1
            pylab.title(multiplot.title)

            # Labels
            pylab.xlabel(multiplot.xlabel, fontsize=fontsize)
            pylab.ylabel(multiplot.ylabel, fontsize=fontsize)

            lines = []
            legends = []
            # List of argument for the plot function
            for j, singleplot in enumerate(multiplot):
                style = singleplot.style
                legend = singleplot.legend
                color = singleplot.color
                label  = singleplot.label
                if(not color):
                    color = self.colors[j % len(self.colors)]

                x = []
                y = []
                labels = []
                for pt in singleplot:
                    x.append(pt.x)
                    y.append(pt.y)

                # no data available. check if label is on.
                if len(x) == 0:
                    if label == True:
                        for i in range(0, singleplot.get_label_size()):
                            x = singleplot.get_label_x(i)
                            y = singleplot.get_label_y(i)
                            labels = singleplot.get_label_text(i)
                            matplotlib.pyplot.text(x, y, labels)
                            pylab.hold(True)
                #        break # nothing else to be done in principle
                    else:
                        print "Warning. Empty data."
                        return
                else:
                    pass
                    # continue to the normal plots


                    # Manage style
                    pointstyle = self.pointstyles[j % len(self.pointstyles)]

                    if "impulses" in style:
                        l = pylab.vlines(x, 0, y, **line2d)
                    elif "linespoints" in style:
                        #l = pylab.plot(x, y, '-', x, y, pointstyle)
                        l = pylab.plot(x, y, pointstyle + '-',**line2d)
                    elif "points" in style:
                        l = pylab.plot(x, y, pointstyle,**line2d)
                    elif "lines" in style:
                        l = pylab.plot(x, y, '-', **line2d)
                    elif "histeps" in style:
                        l = pylab.plot(x, y,  linestyle='steps-mid', **line2d)
                    else:
                        l = pylab.plot(x, y, style,**line2d)

                    if(color):
                        pylab.setp(l, color=color)

                    lines.append(l)
                    legends.append(legend)
            # Legend
            try:
                import warnings
                with warnings.catch_warnings(record=True) as w:
                    #Plot(seq1, Viewpoint="Data", nbcol=2, 
                    #    legend_kwds={'prop':{'size':9}})
                    kwds = {}
                    kwds['prop'] = {'size':legend_size}
                    kwds['ncol'] = legend_nbcol
                    kwds['loc'] = legend_loc
                    if legend_on is True and len(legends)<15:
                        lg = pylab.legend(lines, legends, **kwds)
                        lg.legendPatch.set_alpha(0.1) # transparency
                        if w:
                            # it seems that matplotlib.collections.LineCollection appear
                            # last in legends
                            new_legends = [] # permutation of legends with LineCollection at the end
                            index_begin_line_collection = 0 # where LineCollection objects begin 
                            for li in range(len(lines)):
                                current_legend = legends[li]
                                if str(lines[li].__class__) == "<class 'matplotlib.collections.LineCollection'>":
                                    new_legends += [current_legend]
                                else:
                                    new_legends.insert(index_begin_line_collection, current_legend)
                                    index_begin_line_collection += 1
                            pylab.legend(new_legends, **kwds)
            except Exception, e:
                import warnings
                warnings.warn('legend failed:'+str(e))

            # Grid
            pylab.grid(bool(multiplot.grid))

            # Range
            _xrange = multiplot.xrange
            _yrange = multiplot.yrange
            a = pylab.gca()
            if(_xrange.min != _xrange.max):
                a.set_xlim([round(_xrange.min, 5), round(_xrange.max, 5)])

            if(_yrange.min != _yrange.max):
                a.set_ylim([round(_yrange.min, 5), 
                            round(_yrange.max*options_y_maxrange_ratio, 5)])

            # Tics
            xt = round(multiplot.xtics, 5)
            if(xt > 0):
                xmin, xmax = pylab.xlim()
                pylab.xticks(pylab.arange(xmin, xmax, xt))

            yt = round(multiplot.ytics, 5)
            if(yt > 0):
                ymin, ymax = pylab.ylim()
                pylab.yticks(pylab.arange(ymin, ymax, yt))

        if show == True:
            pylab.show()

        return f1

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
