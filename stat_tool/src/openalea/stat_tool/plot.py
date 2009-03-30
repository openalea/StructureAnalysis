""" Plot functions """
__revision__ = "$Id$"

import sys
import os
import glob

DISABLE_PLOT = False


# !!! Do not plot in nosetests !!!
# buildbot cannot close the windwos popped up by the method/function Plot/plot
# So, we test if the command "python setup.py nosetests" has been used.
# Still, using nosetests executable, windows should pop up.
if("nosetests" in sys.argv):
    DISABLE_PLOT = True



class plotter(object):
    """ Abstract base class for all plotter """

    def plot(self, obj, title, groups=[], *args, **kargs):
        """ Plot obj with title """
        
        raise NotImplementedError()


class fakeplot(plotter):

    def plot(self, obj, title, groups=[], *args, **kargs):
        """ Plot obj with title """
        return

    
class gnuplot(plotter):
    """ GNUPlot implementation """

    def __init__(self):
        """ Initialize GnuPlot """
        import Gnuplot
        self.session = Gnuplot.Gnuplot()

    def plot(self, plotable, title, groups=[], *args, **kargs):
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
        nbx = len(multiset)
        
        # For each subplot
        for i, multiplot in enumerate(multiset):
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
                color = singleplot.color

                x = []
                y = []
                for pt in singleplot:
                    x.append(pt.x)
                    y.append(pt.y)
                
                p = Gnuplot.Data(x, y)
                if style:
                    "todo: check that this option works."
                    p.set_option(with_=style)

                if legend:
                    p.set_option(title=legend)
                plot_list.append(p)
                

                
            # Range
            xrange = multiplot.xrange
            yrange = multiplot.yrange
            if(xrange.min != xrange.max):
                g('set xrange[%f:%f]'%(xrange.min, xrange.max))
            if(yrange.min != yrange.max):
                g('set yrange[%f:%f]'%(yrange.min, yrange.max))

            # Tics
            if(multiplot.xtics > 0):
                g('set xtics 0, %f'%(multiplot.xtics))
            if(multiplot.ytics > 0):
                g('set ytics 0, %f'%(multiplot.ytics))


            g.plot(*plot_list)
            raw_input("Press Enter to continue")


class mplotlib(plotter):
    """matplotlib implementation"""

    linestyles = ('-', '--', ':', '.')
    pointstyles = ('o', '^', 'x', '+', 's', 'v', '>', '<')
    colors = ('g', 'r', 'y', 'b', 'k', 'm', 'c')
    
    def __init__(self):
        """ Initialize matplotlib """
        
        import matplotlib
        #matplotlib.use('Qt4Agg')
        import pylab
        self.pylab = pylab

    def plot(self, plotable, title, groups=[], show=True, *args, **kargs):
        """ 
        Plot a plotable with title 
        groups : list of group (int) to plot
        show=True by default will pop up the figure
        """
        
        pylab = self.pylab
        multiset = plotable

        # Title & border
        #if title: pylab.suptitle(title) 

        #multiset.border

        # Configure figure
        f1 = pylab.figure(1)
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

                f = pylab.figure(g+1)
                f.set_facecolor("w")

        # nb subplot
        nbx = len(multiset)
        
        # For each subplot
        for i, multiplot in enumerate(multiset):
            
            g = multiplot.group
            # Group filter
            if(groups and g not in groups):
                continue

            # Select window
            pylab.figure(g+1)

            # Select Subplot
            pylab.subplot(group_size[g], 1, group_index[g] + 1)
            group_index[g] += 1

            pylab.title(multiplot.title)
            
            # Labels
            pylab.xlabel(multiplot.xlabel)
            pylab.ylabel(multiplot.ylabel)
            
            lines = []
            legends = []
            # List of argument for the plot function
            for j, singleplot in enumerate(multiplot):

                style = singleplot.style 
                legend = singleplot.legend
                color = singleplot.color

                if(not color): color = self.colors[j%len(self.colors)]

                x = []
                y = []
                for pt in singleplot:
                    x.append(pt.x)
                    y.append(pt.y)
                
                # Manage style
                pointstyle = self.pointstyles[j%len(self.pointstyles)]

                if "impulses" in style:
                    l = pylab.vlines(x, 0, y)
                elif "lines" and  "points" in style:
                    #l = pylab.plot(x, y, '-', x, y, pointstyle)
                    l = pylab.plot(x, y, pointstyle + '-')
                elif "points" in style:
                    l = pylab.plot(x, y, pointstyle)
                elif "lines" in style:
                    l = pylab.plot(x, y, '-')
                else:
                    l = pylab.plot(x, y, style)
                
                if(color):
                    pylab.setp(l, color=color)
                    
                lines.append(l)
                legends.append(legend)

            # Legend
            lg = pylab.legend(lines, legends)
            lg.legendPatch.set_alpha(0.1) # transparency

            # Grid
            pylab.grid(bool(multiplot.grid))
                
            # Range
            xrange = multiplot.xrange
            yrange = multiplot.yrange
            a = pylab.gca()
            if(xrange.min != xrange.max):
                a.set_xlim([round(xrange.min, 5), round(xrange.max, 5)])

            if(yrange.min != yrange.max):
                a.set_ylim([round(yrange.min, 5), round(yrange.max,5)])

            # Tics
            xt = round(multiplot.xtics, 5)
            if(xt > 0):
                xmin, xmax = pylab.xlim() 
                pylab.xticks(pylab.arange(xmin, xmax, xt))

            yt = round(multiplot.ytics, 5)
            if(yt > 0):
                ymin, ymax = pylab.ylim() 
                pylab.yticks(pylab.arange(ymin, ymax, yt))

        if show==True:
            pylab.show()
   

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
    if(PLOTTER):
        return PLOTTER

    # Try to import a plotter
    else:
        try:
            plotter = mplotlib()
            
        except ImportError:
            plotter = gnuplot()

        set_plotter(plotter)
        return plotter
