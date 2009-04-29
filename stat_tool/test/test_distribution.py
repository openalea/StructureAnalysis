"""Distribution tests"""
__revision__ = "$Id$"


from openalea.stat_tool import _stat_tool
from openalea.stat_tool.distribution import Distribution, Uniform, Binomial
from openalea.stat_tool.distribution import NegativeBinomial, Poisson
from openalea.stat_tool.distribution import ToHistogram, ToDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.output import Display,  Plot
from openalea.stat_tool import Estimate
from openalea.stat_tool import Simulate

from openalea.stat_tool import Cluster, \
    Transcode, ValueSelect, Shift, Compare, ComparisonTest, Fit


from tools import interface

class Test(interface):
    """a simple unittest class

 
    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/distribution1.dist",
                           Distribution)
    
    def build_data(self):
        d1 = Binomial(0, 10, 0.5)
        d2 = Distribution("BINOMIAL", 0, 10, 0.5)
        assert d1 == d2
        return d1

    def test_empty(self):
        self.empty()

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()
        
    def test_constructors(self):
        h = Histogram([1,2,3,4,5,6,1,2,3])
        assert h
        
        # from histogram
        dist = Distribution(h)
        assert dist
        
        #from parametric model
        pm = _stat_tool._ParametricModel(h)
        dist = Distribution(pm)
        assert dist
        

    def test_print(self):
        self.print_data()
        
    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()
        
    def test_len(self):
        """not implemented; irrelevant?"""
        pass
    
    def test_plot(self):        
        self.plot()

    def test_save(self):
        self.save()
                
    def test_plot_write(self):
        self.plot_write()

    def test_file_ascii_write(self):
        self.file_ascii_write()
      
    def test_spreadsheet_write(self):
        self.spreadsheet_write()
    
    def test_simulate(self):
        self.simulate()
        
    def test_survival_ascii_write(self):
        self.survival_ascii_write()
        
    def test_survival_spreadsheet_write(self):
        self.survival_spreadsheet_write()
        
    def test_extract(self):
        pass

    def test_extract_data(self):
        s = self.simulate()
        e = s.estimate_parametric("B")
        d = e.extract_data()
        assert d
        _eprime = Estimate(s, "Binomial")
        
    def test_truncate(self):
        # todo:find a test
        s = self.data
        _res = s.truncate(4) 
        

class TestDistribution:
    """Test the distribution (Uniform, Binomial, ...)
    
    test the sup_bound, inf_bound, probability, parameter,ident
    
    test the ToDistribution and ToHistogram
    """
    def test_to_histogram(self):

        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        h = d.simulate(1000)
        d2 = ToDistribution(h)
        assert h and d2

        h2 = ToHistogram(d2)
        assert h2
        assert h == h2

    def test_uniform(self):

        d = Distribution("UNIFORM", 0, 10)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == -1
        assert d.get_ident == 4
        
        d = Uniform(0, 10)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)

    def test_binomial(self):

        d = Distribution("BINOMIAL", 0, 10, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == 10
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == -1
        assert d.get_ident == 1

        d = Binomial(0, 10, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)

    def test_poisson(self):

        d = Distribution("POISSON", 0, 2)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == -1
        assert d.get_parameter == 2
        assert d.get_ident == 2

        d = Poisson(0, 2)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)

    def test_neg_binomial(self):
        
        d = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.5)
        assert list(d.simulate(1000))
        assert d.get_sup_bound == -1
        assert d.get_inf_bound == 0
        assert d.get_probability == 0.5
        assert d.get_parameter == 1
        assert d.get_ident == 3
        d = NegativeBinomial(0, 1, 0.5)
        assert list(d.simulate(1000))

        m = d.simulate(1000).extract_model()
        assert isinstance(m, _stat_tool._ParametricModel)
        
    def test_multimodial(self):
        pass
    
    def test_getters(self):
        dist = Distribution(Histogram([1, 1, 1, 2, 2, 2, 3, 3, 3]))
        assert dist.get_mean == 2
        assert 0.3333 < dist.get_max < 0.3334 
        assert dist.get_mean == 2
        assert 0.6666 < dist.get_variance < 0.6667


def test1():
    #########################################################################
    #
    #  discrete distributions/histograms, comparison of histograms
    # /frequency distributions
    #
    #  beech, Wild cherry tree: number of nodes per growth unit (GU)
    #
    #  meri1.his: order 1,
    #  meri1.his: order 2,
    #  meri1.his: order 3, GU 1,
    #  meri1.his: order 3, GU 2,
    #  meri5.his: short shoots.
    #
    #########################################################################


    _dist0 = Distribution("NEGATIVE_BINOMIAL", 0, 1, 0.3)
    dist0 = Distribution("data/distribution1.dist")

    dist1 = Distribution("B", 0, 10, 0.3)
    dist1 = Distribution("NB", 0, 3.5, 0.3)

    # computation of survivior function and hasard rate

    Display(dist1, ViewPoint="Survival")
    Plot(dist1, ViewPoint="Survival")

    # simulation/estimation

    histo1 = Simulate(dist1, 200)

    #Display(histo1, ViewPoint->"Survival")
    #Plot(histo1, ViewPoint->"Survival")

    # "B" "BINOMIAL", "P" "POISSON", "NB" "NEGATIVE_BINOMIAL"
    # InfBoundStatus->"Free" (default) / "Fixed"
    # MinInfBound->0 (default) / 1

    _dist2 = Estimate(histo1, "NB", MinInfBound=0, InfBoundStatus="Fixed")

    fagus = Histogram("data/fagus1.his")

    # transformation of histograms, extraction/filter

    _histo2 = Cluster(fagus, "Step", 2)
    _histo3 = Cluster(fagus, "Information", 0.8)
    _histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])
    histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
    Display(histo5, Detail=2)

    _histo7 = Shift(fagus, -2)

    _histo8 = ValueSelect(fagus, 2, 8)

    _dist3 = Estimate(fagus, "B")

    # comparison of histograms

    meri1 = Histogram("data/meri1.his")
    meri2 = Histogram("data/meri2.his")
    meri3 = Histogram("data/meri3.his")
    meri4 = Histogram("data/meri4.his")
    meri5 = Histogram("data/meri5.his")

    ## Compare(meri1, meri2, meri3, meri4, meri5, "N", FileName="ASCII/meri.cmp")
    Compare(meri1, meri2, meri3, meri4, meri5, "O")

    ComparisonTest("F", meri1, meri2)
    ComparisonTest("T", meri1, meri2)
    ComparisonTest("W", meri1, meri2)

    # fit of a known distribution to an  histogram

    dist5 = Fit(meri5, Distribution("B", 0, 10, 0.437879))
    Display(dist5, Detail=2)
    Plot(dist5)
    



    
    
