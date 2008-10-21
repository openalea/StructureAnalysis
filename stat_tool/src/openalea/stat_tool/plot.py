__docformat__ = "restructuredtext"
__doc__ = """ Plot functions """


import sys

DISABLE_PLOT = False


# !!! Do not plot in nosetests !!!
if("nosetests" in sys.argv[0]):
    DISABLE_PLOT = True



class plotter(object):
    """ Abstract base class for all plotter """

    def plot(self, obj, title, groups=[], *args, **kargs):
        """ Plot obj with title """
        
        raise NotImplementedError()


import os, glob

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
                if(style): p.set_option(with=style)
                if(legend) : p.set_option(title=legend)
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
    """ 
    MathPlotLib implementation 
    """

    linestyles = ('-', '--', ':', '.')
    pointstyles = ('o', '^', 'x', '+', 's', 'v', '>', '<')
    colors = ('g', 'r', 'y', 'b', 'k', 'm', 'c')
    

    def __init__(self):
        """ Initialize matplotlib """
        
        import matplotlib
        #matplotlib.use('Qt4Agg')

        import pylab
        self.pylab = pylab



    def plot(self, plotable, title, groups=[], *args, **kargs):
        """ 
        Plot a plotable with title 
        groups : list of group (int) to plot
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

                  


        pylab.show()


################################################################################    

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




########################## Test ################################################

class Test:
    
    def test_plotable(self):
        import _stat_tool


        p = _stat_tool.SinglePlot()
        p.add_point(_stat_tool.PlotPoint(1.0, 1.2))
        p.add_point(_stat_tool.PlotPoint(2.0, 3.2))

        assert len(list(p)) == 2

        for i in p:
            assert (i.x, i.y)

        assert list(p)

        p = _stat_tool.MultiPlot(2)
        assert p
        assert len(p) == 2


        p = _stat_tool.MultiPlotSet(3)
        assert len(p) == 3
        


    def test_plotable2(self):
        import _stat_tool

        p = _stat_tool.MultiPlotSet(1)
        p[0][0].add_point(_stat_tool.PlotPoint(1.0, 1.2))
        p[0][0].add_point(_stat_tool.PlotPoint(3.0, 5.2))
                
        assert len(p) == 1

        # Test iterator
        multiset = p
        for i, multiplot in enumerate(multiset):
            assert multiplot

            for singleplot in multiplot:
                assert singleplot and len(multiplot)

        # test plot
        plotter = get_plotter()
        plotter.plot(p, "test_plot")
        

        
    def get_plotable(self):
        import _stat_tool

        p = _stat_tool.MultiPlotSet(3)
        p.title = "TestPlot"
                
        p[0].title = "P1"
        p[0].xlabel = "x1"
        p[0].ylabel = "y1"
                
        p[0][0].add_point(_stat_tool.PlotPoint(1.0, 1.2))
        p[0][0].add_point(_stat_tool.PlotPoint(3.0, 5.2))
        p[0][0].style = "lines"


        p[1].resize(2)
                
        p[1].title = "P2"
        p[1].xlabel = "x2"
        p[1].ylabel = "y2"
        
        p[1][0].add_point(_stat_tool.PlotPoint(12.0, 1.2))
        p[1][0].add_point(_stat_tool.PlotPoint(31.0, 5.2))
        p[1][0].legend = "PLOT1"
        p[1][0].color = "g"
        p[1][0].style = "points"
        

        p[1][1].add_point(_stat_tool.PlotPoint(14.0, 5.2))
        p[1][1].add_point(_stat_tool.PlotPoint(5.0, 2))
        p[1][1].add_point(_stat_tool.PlotPoint(3.0, 5))
        p[1][1].add_point(_stat_tool.PlotPoint(35.0, 3.2))
        p[1][1].legend = "PLOT2"
        p[1][1].style = "impulses"
                
        p[2].title = "P3"
        p[2].xlabel = "x3"
        p[2].ylabel = "y3"
        p[2].xtics = 1.5
        p[2].group = 1


        p[2][0].add_point(_stat_tool.PlotPoint(3, 0.2))
        p[2][0].add_point(_stat_tool.PlotPoint(8, 7.2))
        p[2][0].add_point(_stat_tool.PlotPoint(12.0, 1.2))
        p[2][0].add_point(_stat_tool.PlotPoint(31.0, 5.2))
        p[2][0].legend = "PLOT3"
        p[2][0].color = "y"
        p[2][0].style = "linespoints"
        


        assert len(p) == 3

        return p

    def test_matplotlib(self):

        a = self.get_plotable()
        plotter = get_plotter()
        plotter.plot(a, "test_plot")


    def test_gnuplot(self):
        
        if(DISABLE_PLOT) : return
        a = self.get_plotable()
        plotter = gnuplot()
        plotter.plot(a, "test_plot")

        
