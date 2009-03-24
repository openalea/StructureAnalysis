"""output tests"""
__revision__ = "$Id: $"


from openalea.stat_tool.distribution import Uniform, Distribution
from openalea.stat_tool.histogram import Histogram 
from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.convolution import Convolution
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool.data_transform import Shift
from openalea.stat_tool.output import Plot


class Test:
    """a simple unittest class"""
    
    def get_mixture(self):
        """create a mixture data""" 
        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)
        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        return m

    def get_mixture_2(self):
        """create another mixture data""" 
        h = Histogram(("meri2.his"))
        m = h.estimate_mixture(["B", "NB"])
        return m

    def test_old_plot(self):
        m = self.get_mixture()
        m.old_plot()

    def test_plot_mixture_1(self):
        m = self.get_mixture()
        m.plot()

    def test_plot_mixture_2(self):
        m = self.get_mixture_2()
        m.plot()

    def test_plot_mixture_data(self):
        mixt1 = Mixture(0.6, Distribution("B", 2, 18, 0.5), 0.4, 
                        Distribution("NB", 10, 10, 0.5))
        mixt_histo1 = Simulate(mixt1, 200)

        mixt1.plot()
        mixt_histo1.plot()

    def test_plot_convolution(self):
        
        convol1 = Convolution(("convolution1.conv"))
        convol1.plot()

        histo_b2 = Histogram(("nothofagus_antarctica_bud_2.his"))
        histo_s2 = Histogram(("nothofagus_antarctica_shoot_2.his"))
        
        convol31 = Estimate(Shift(histo_s2, 1), "CONVOLUTION", Estimate(histo_b2, "NP"), 
                            NbIteration=100, Estimator="PenalizedLikelihood", Weight=0.5)
        convol31.plot()

    def test_plot_convolution_data(self):
        
        convol1 = Convolution(("convolution1.conv"))
        convol_histo1 = Simulate(convol1, 200)
        convol_histo1.plot()

    def test_plot_distribution_set(self):
        
        d1 = Distribution("B", 2, 18, 0.5) 
        d2 = Distribution("NB", 10, 10, 0.5)
        d3 = Distribution("U", 10, 20)
        
        Plot(d1, d2, d3)
        d1.old_plot()
        
    def test_plot_survival(self):
        

        d1 = Distribution("B", 2, 18, 0.5) 
        
        d1.plot(ViewPoint="survival")

        histo1 = Simulate(d1, 200)
        histo1.plot(ViewPoint="survival")
        
    def test_plot_parametric_model(self):
        
        dist1 = Distribution("NB", 0, 3.5, 0.3)
        histo1 = Simulate(dist1, 200)
        Plot(histo1)
        dist2 = Estimate(histo1, "NB", MinInfBound=0, InfBoundStatus="Fixed")
        Plot(dist2) 

def test1():
    pass    



if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    
    test1()
