"""output tests"""
__version__ = "$Id: test_output.py 9397 2010-08-10 10:50:07Z cokelaer $"

from openalea.stat_tool.distribution import Uniform, Distribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.data_transform import Shift
from openalea.stat_tool.output import Plot, Display
from openalea.stat_tool.plot import DISABLE_PLOT

from tools import runTestClass

from openalea.stat_tool import get_shared_data

class Test:
    """a simple unittest class"""

    def __init__(self):
        pass

    def get_mixture(self):
        """create a mixture data"""
        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)
        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        return m

    def get_mixture_2(self):
        """create another mixture data"""
        h = Histogram(get_shared_data("meri2.his"))
        m = h.estimate_mixture("B", "NB")
        return m

    def _test_old_plot(self):
        m = self.get_mixture()
        if DISABLE_PLOT == False:
            m.old_plot()

    def test_plot_mixture_1(self):
        m = self.get_mixture()
        if DISABLE_PLOT == False:
            m.plot()

    def test_plot_mixture_2(self):
        m = self.get_mixture_2()
        if DISABLE_PLOT == False:
            m.plot()

    def test_plot_mixture_data(self):
        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 0.4,
                        Distribution("NB", 10, 10, 0.5))
        mixt_histo1 = Simulate(mixt1, 200)

        if DISABLE_PLOT == False:
            mixt1.plot()
            mixt_histo1.plot()

    def test_plot_convolution(self):

        convol1 = Convolution("data/convolution1.conv")
        if DISABLE_PLOT == False:
            Plot(convol1.extract_elementary(1), convol1.extract_elementary(2))

        histo_b2 = Histogram("data/nothofagus_antarctica_bud_2.his")
        histo_s2 = Histogram("data/nothofagus_antarctica_shoot_2.his")

        convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION",
                            Estimate(histo_b2, "NP"),
                            NbIteration=100,
                            Estimator="PenalizedLikelihood",
                            Weight=0.5)
        if DISABLE_PLOT == False:
            Plot(convol31.extract_elementary(1))

    def test_plot_convolution_data(self):

        convol1 = Convolution("data/convolution1.conv")
        convol_histo1 = Simulate(convol1, 200)
        if DISABLE_PLOT == False:
            convol_histo1.plot()

    def _test_plot_distribution_set(self):

        d1 = Distribution("B", 2, 18, 0.5)
        d2 = Distribution("NB", 10, 10, 0.5)
        d3 = Distribution("U", 10, 20)

        if DISABLE_PLOT == False:
            Plot(d1, d2, d3)
            d1.old_plot()

    def test_plot_survival(self):


        d1 = Distribution("B", 2, 18, 0.5)

        if DISABLE_PLOT == False:
            d1.plot(ViewPoint="Survival")

        histo1 = Simulate(d1, 200)
        if DISABLE_PLOT == False:
            histo1.plot(ViewPoint="Survival")

    def test_plot_parametric_model(self):

        dist1 = Distribution("NB", 0, 3.5, 0.3)
        histo1 = Simulate(dist1, 200)
        if DISABLE_PLOT == False:
            Plot(histo1)
        dist2 = Estimate(histo1, "NB", MinInfBound=0, InfBoundStatus="Fixed")
        if DISABLE_PLOT == False:
            Plot(dist2)



def test_output_plot_viewpoint_survival():
    d1 = Distribution("B", 2,19,0.5)
    Plot(d1)
    Plot(d1, ViewPoint='Survival')

def test_output_display_detail():
    d1 = Distribution("B", 2,19,0.5)
    print Display(d1)
    print Display(d1, Detail=1)
    print Display(d1, Detail=2)
    try:
        print Display(d1, Detail=3)
        assert False
    except:
       assert ValueError

def test_output_display_viewpoint_survival():
    d1 = Distribution("B", 2,19,0.5)
    print Display(d1, ViewPoint="Survival")
    convol1 = Convolution("data/convolution1.conv")
    print Display(convol1, ViewPoint="Survival")

def test_output_display_viewpoint_data():
    from openalea.stat_tool import Vectors
    v = Vectors('data/vectors.vec')
    Display(v, ViewPoint='Data') 
    Display(v) 
    #ddata viewpoint is dedicated to Vectors or sequences
    try:
        d1 = Distribution("B", 2,19,0.5)
        Display(d1, ViewPoint="Data")
        assert False
    except:
        assert True

if __name__ == "__main__":
    runTestClass(Test())

