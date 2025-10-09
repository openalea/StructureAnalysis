"""output tests"""

__version__ = "$Id$"

from openalea.stat_tool import _stat_tool
from openalea.stat_tool.plot import DISABLE_PLOT, get_plotter, gnuplot

DISABLE_PLOT = True

from .tools import runTestClass


class Test:
    def __init__(self):
        pass

    def test_plotable(self):
        p = _stat_tool.SinglePlot()
        p.add_point(_stat_tool.PlotPoint(1.0, 1.2))
        p.add_point(_stat_tool.PlotPoint(2.0, 3.2))

        import openalea

        assert isinstance(p, openalea.stat_tool._stat_tool.SinglePlot)

        p2 = _stat_tool.MultiPlot(2)
        assert p2
        assert len(p2) == 2

        p = _stat_tool.MultiPlotSet(3)
        assert len(p) == 3

    def test_plotable2(self):
        # TODO: understand why the three points do not show up.
        p = _stat_tool.MultiPlotSet(1)

        p[0].resize(1)
        p[0].set_xrange((-10, 10))
        p[0].set_yrange((-5, 8))
        p[0].xlabel = "x"
        p[0].ylabel = "y"
        # equivalent to p[0].xrange = _stat_tool.PlotPoint((0 , 3.5)
        p[0].xtics = 1
        p[0].ytics = 1

        p[0][0].add_point(_stat_tool.PlotPoint(1.0, 1.2))
        p[0][0].add_point(_stat_tool.PlotPoint(3.0, 5.2))
        p[0][0].add_point(_stat_tool.PlotPoint(2.0, 2.3))
        p[0][0].add_text(0.5, 0.2, "*")
        p[0][0].add_text(3.5, 5.5, "*")
        p[0][0].legend = "legend"
        p[0][0].color = "r"
        p[0][0].style = "+"
        # p[0][0].style = "linespoints";
        # p[0][0].label = True

        from openalea.stat_tool.output import plot

        plotter = plot.get_plotter()
        plotter.plot(p, "Points")

        assert len(p) == 1

        import openalea

        # Test iterator
        multiset = p
        for _i, multiplot in enumerate(multiset):
            assert multiplot

            for singleplot in multiplot:
                assert isinstance(
                    singleplot, openalea.stat_tool._stat_tool.SinglePlot
                ) and len(multiplot)

        # test plot
        plotter = get_plotter()
        if DISABLE_PLOT == False:
            plotter.plot(p, "test_plot")

    def get_plotable(self):
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
        if DISABLE_PLOT == False:
            plotter.plot(a, "test_plot")

    def test_gnuplot(self):
        if DISABLE_PLOT:
            return
        a = self.get_plotable()
        plotter = gnuplot()
        if DISABLE_PLOT == False:
            plotter.plot(a, "test_plot")


if __name__ == "__main__":
    runTestClass(Test())
