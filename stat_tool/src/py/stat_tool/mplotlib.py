from . import __stat_tool
from matplotlib import pyplot as plt
from functools import wraps
import numpy as np
from types import ModuleType

linestyles = ('-', '--', ':', '.')
pointstyles = ('o', '^', 'x', '+', 's', 'v', '>', '<')
colors = ('g', 'r', 'y', 'b', 'k', 'm', 'c')

def wrapper(f):
    @wraps(f)
    def plot(self, title=None, groups=None, *args, **kargs):
        """Plot a plotable with title

        :param plotable: a plotable instance from standard stat_tool objects such as :func:`Distribution`
        :param groups: list of group (int) to plot
        :param nbcol:
        :param legend_size: 10
        :param legend_nbcol: 2
        :param legend_loc: best
        :param legend: True/False

        y_maxrange_ratio=1 multiply max y range by this value
        """

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

        multiset = f(self)

        if len(multiset)==1:
            nbcol = 1
        # Title & border
        #if title: pylab.suptitle(title)
        #multiset.border

        # Configure figure
        fig = plt.figure(fig_id, figsize=(10, 10))
        #fig = plt.gcf()
        fig.clf()

        fig.set_facecolor("w")

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

                _fig = plt.figure(fig_id + g )
                _fig.set_facecolor("w")

        # nb subplot
        _nbx = len(multiset)

        if title:
            fig.suptitle(title, fontsize=fontsize)
        # For each subplot
        for i, multiplot in enumerate(multiset):
            g = multiplot.group
            # Group filter
            if(groups and g not in groups):
                continue

            # Select window
            plt.figure(fig_id + g )

            # Select Subplot
            axes = plt.subplot(group_size[g]%nbcol+group_size[g]/nbcol,
                          nbcol , group_index[g] + 1)
            group_index[g] += 1
            axes.set_title(multiplot.title)

            # Labels
            axes.set_xlabel(multiplot.xlabel, fontsize=fontsize)
            axes.set_ylabel(multiplot.ylabel, fontsize=fontsize)

            lines = []
            legends = []
            # List of argument for the plot function
            for j, singleplot in enumerate(multiplot):
                style = singleplot.style
                legend = singleplot.legend
                color = singleplot.color
                label  = singleplot.label
                if(not color):
                    color = colors[j % len(colors)]

                x = []
                y = []
                labels = []
                for pt in singleplot:
                    x.append(pt.first)
                    y.append(pt.second)

                # no data available. check if label is on.
                if len(x) == 0:
                    if label == True:
                        for i in range(0, len(singleplot)):
                            x = singleplot[i].x
                            y = singleplot[i].y
                            labels = singleplot[i].label
                            plt.text(x, y, labels)
                            fig.hold(True)
                #        break # nothing else to be done in principle
                    else:
                        print("Warning. Empty data.")
                        return
                else:
                    pass
                    # continue to the normal plots


                    # Manage style
                    pointstyle = pointstyles[j % len(pointstyles)]

                    if "impulses" in style:
                        l = plt.vlines(x, 0, y, color=color, **line2d)
                    elif "linespoints" in style:
                        #l = pylab.plot(x, y, '-', x, y, pointstyle)
                        l = plt.plot(x, y, pointstyle + '-', color=color,**line2d)
                    elif "points" in style:
                        l = plt.plot(x, y, pointstyle, color=color,**line2d)
                    elif "lines" in style:
                        l = plt.plot(x, y, '-', color=color, **line2d)
                    elif "histeps" in style:
                        l = plt.plot(x, y, color=color, linestyle='steps-mid', **line2d)
                    else:
                        l = plt.plot(x, y, style, color=color, **line2d)

                    #if(color):
                    #    l.set_color(color)

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
                        lg = axes.legend(lines, legends, **kwds)
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
                            axes.legend(new_legends, **kwds)
            except Exception as e:
                import warnings
                warnings.warn('legend failed:'+str(e))

            # Grid
            axes.grid(bool(multiplot.grid))

            # Range
            _xrange = multiplot.xrange
            _yrange = multiplot.yrange
            #a = pylab.gca()
            if(_xrange.first != _xrange.second):
                axes.set_xlim(sorted([round(_xrange.first, 5), round(_xrange.second, 5)]))

            if(_yrange.first != _yrange.second):
                ymin, ymax = sorted([_yrange.first, _yrange.second])
                ymin = round(ymin, 5)
                ymax = round(ymax * options_y_maxrange_ratio, 5)
                axes.set_ylim([ymin, ymax])

            # Tics
            xt = round(multiplot.xtics, 5)
            if(xt > 0):
                xmin, xmax = axes.get_xlim()
                axes.set_xticks(np.arange(xmin, xmax, xt))

            yt = round(multiplot.ytics, 5)
            if(yt > 0):
                ymin, ymax = axes.get_ylim()
                axes.set_yticks(np.arange(ymin, ymax, yt))

        return fig

    return plot

def get_plotables(obj):
    if isinstance(obj, ModuleType):
        for obj in obj.__dict__.itervalues():
            get_plotables(obj)
    elif hasattr(obj, 'get_plotable'):
        obj.plot = wrapper(getattr(obj, 'get_plotable'))
        #delattr(obj, 'get_plotable')

def load_ipython_extension(ipython):

    get_plotables(__stat_tool.stat_tool)

def unload_ipython_extension(ipython):

    pass
