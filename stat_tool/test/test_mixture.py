"""Mixture tests
#
#  Frequency distributions
#
#  Objective: Analyzing the number of nodes of growth units in selected architectural
#             position considering the respective roles of preformation and neoformation,
#
#  Methods: comparison tests, one-way variance analysis,
#           estimation of finite mixture of distributions.
#
#  Wild cherry tree: number of nodes per growth unit (GU)
#
#  Data: Dominique Fournier
#
#  meri1.his: order 1,
#  meri1.his: order 2,
#  meri1.his: order 3, GU 1,
#  meri1.his: order 3, GU 2,
#  meri5.his: short shoots.
#
#
#  Poplar: number of nodes per growth unit
#
#  Data: Yves Caraglio and Herve Rey
#
#  peup1.his: order 2,
#  peup2.his: order 3,
#  peup3.his: order 4,
#  peup4.his: order 5,
#  peup5.his: order 3, GU 4,
#  peup6.his: order 3, acrotony.
#
#########################################################################
"""
__revision__ = "$Id: $"

import os


from openalea.stat_tool.distribution import Uniform
from openalea.stat_tool.mixture import Mixture
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.mixture import _MvMixture
from openalea.stat_tool.distribution import Binomial, Poisson
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.comparison import Compare, ComparisonTest

class Test:
    """a simple unittest class"""

    def test_empty(self):
        """Test that empty constructor fails"""
        try:
            m = Mixture()
            assert False
        except Exception:
            assert True

    def test_file(self):
        """run constructor with filename argument"""
        m = Mixture("mixture1.mixt")
        assert m

    def test_build_mixture(self):
        """run constructor with two distributions as arguments"""

        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        
        assert m
        return m

    def test_plot(self):
        """run plotting routines """
        m = self.test_build_mixture()
        m.plot()

        # todo move to matplotlib
        assert str(m)
        m.display()

    def test_simulation(self):
        """Test the simulate method"""
        m = self.test_build_mixture()
        s = m.simulate(10)

        assert s.nb_component() == 3
        assert str(s)

    def test_extract(self):
        """run and test the extract methods"""

        m = self.test_build_mixture()
        assert m.nb_component() == 3

        assert m.extract_weight() == ExtractDistribution(m, "Weight")

        assert m.extract_mixture() == ExtractDistribution(m, "Mixture")

        assert ExtractDistribution(m, "Component", 1) == Uniform(0, 10)
        assert ExtractDistribution(m, "Component", 2) == Uniform(10, 20)
        assert ExtractDistribution(m, "Component", 3) == Uniform(20, 30)

        assert m.extract_component(1) == Uniform(0, 10)
        assert m.extract_component(2) == Uniform(10, 20)
        assert m.extract_component(3) == Uniform(20, 30)

    def test_extract_data(self):
        """run and test the extract_data methods"""

        h = Histogram("meri2.his")
        m = h.estimate_mixture(["B", "NB"])

        d = m.extract_data()
        assert d

    def test_build_mv_mixture(self):

        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(0, 12, 0.5)
        d13 = Binomial(0, 12, 0.8)

        d21 = Poisson(0, 18.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, .20)

        m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        assert m
        return m

    def _test_mv_fromfile(self):
        # From file
        m = _MvMixture("mixture_mv1.mixt")
        assert m

        # File
        m.save("test_mv.mixt")

        _m1 = _MvMixture("mixture_mv1.mixt")
        m2 = _MvMixture("test_mv.mixt")
        assert m.nb_component() == m2.nb_component()
        assert str(m) == str(m2)

        os.remove("test_mv.mixt")

        mnp = _MvMixture("mixture_mv_nonparam.mixt")
        assert mnp

        try:
            _h = _MvMixture("no_such_file.mixt")
            assert False
        except Exception:
            assert True

    def test_simulate_estimate_mv_mixture(self):

        d11 = Binomial(0, 12, 0.1)
        d12 = Binomial(2, 13, 0.6)
        d13 = Binomial(3, 15, 0.9)

        d21 = Poisson(0, 25.0)
        d22 = Poisson(0, 5.0)
        d23 = Poisson(0, 0.2)


        m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
        v = m.simulate(5000)
        assert v

        m_estim_model = v.mixture_estimation(m, 100, [True, True])
        assert m_estim_model
        m_estim_nbcomp = v.mixture_estimation(2)
        assert m_estim_nbcomp

#  funtional tests
from openalea.stat_tool.output import plot, Plot
from openalea.stat_tool.data_transform import Merge, Shift, ExtractData
from openalea.stat_tool.simulate import Simulate
import openalea.stat_tool.distribution as distribution 
from openalea.stat_tool.cluster import Cluster

def _test1():
    plot.DISABLE_PLOT = True
 
    meri1 = Histogram("meri1.his")
    meri2 = Histogram("meri2.his")
    meri3 = Histogram("meri3.his")
    meri4 = Histogram("meri4.his")
    meri5 = Histogram("meri5.his")

    Plot(meri1, meri2, meri3, meri4, meri5)
    Compare(meri1, meri2, meri3, meri4, meri5, "N")


    ComparisonTest("F", meri1, meri2)
    ComparisonTest("T", meri1, meri2)
    ComparisonTest("W", meri1, meri2)

    ComparisonTest("F", meri1, meri3)
    ComparisonTest("T", meri1, meri3)
    ComparisonTest("W", meri1, meri3)

    # estimation of a mixture of two distributions assuming a first sub-population of GUs
    # made only of a preformed part and a second sub-population made of both a preformed part
    # and a neoformed part
    
    _mixt1 = Estimate(meri2, "MIXTURE", "B", "B")
    
    meri = Merge(meri1, meri2, meri3, meri4, meri5)
    
    # model selection approach: estimation of both the mixture parameters and
    # the number of components 
    
    mixt2 = Estimate(meri, "MIXTURE", "B", "B", "B", "B",  NbComponent="Estimated")
    # mixt2 = Estimate(meri, "MIXTURE", "NB", "NB")
    # Plot(ExtractDistribution(mixt2, "Mixture"))
    # Display(mixt2)
    
    _mixt_data = ExtractData(mixt2)
    
    
    dist5 = Estimate(meri5, "BINOMIAL")
    # Display(dist5, Detail->2)
    # Plot(dist5)
    
    histo5 = Simulate(dist5, 100)
    # Display(histo5, Detail->2)
    # Plot(histo5)
    
    
    peup1 = Histogram("peup1.his")
    peup2 = Histogram("peup2.his")
    peup3 = Histogram("peup3.his")
    peup4 = Histogram("peup4.his")
    peup5 = Histogram("peup5.his")
    peup6 = Histogram("peup6.his")
    
    mixt10 = Estimate(peup2, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    
    peup = Merge(peup1, peup2, peup3, peup4, peup5, peup6)
    
    histo1 = Shift(peup, -1)
    histo2 = Cluster(peup, "Information", 0.8)
    histo3 = Cluster(peup, "Step", 10)
    histo4 = Cluster(peup, "Limit", [13, 24])
    # Display(histo4, Detail->2)
    # Plot(histo4)


    mixt11 = Estimate(peup, "MIXTURE", "B", "NB", "NB", "NB", NbComponent="Estimated")
    # mixt11 = Estimate(peup, "MIXTURE", "B", "NB")

    d11 = distribution.Binomial(0, 12, 0.1)
    d12 = distribution.Binomial(2, 13, 0.6)
    d13 = distribution.Binomial(3, 15, 0.9)
    
    d21 = distribution.Poisson(0, 25.0)
    d22 = distribution.Poisson(0, 5.0)
    d23 = distribution.Poisson(0, 0.2)

    m = _MvMixture([0.1, 0.2, 0.7], [[d11, d21], [d12, d22], [d13, d23]])
    print m

    m2 = _MvMixture("mixture_mv1.mixt")
    print m2

    print "Egalite des melanges construits par liste ",\
      "de distributions et par fichiers : ", str(str(m)==str(m2))

    m = _MvMixture("mixture_mv_nonparam.mixt")
    print m

    print "Simulation de melanges multivaries : "
    v = m.simulate(5000)
    print v

    m.plot(variable=1, Title="Simulated mixture")

    print "Estimation de melanges multivaries ", \
        "d'apres un modele initial : "
    m_estim_model = v.mixture_estimation(m, 100,  [True, True])
      
    m_estim_model.plot(variable = 1, Title="Estimated mixture")
    
    print "Estimation de melanges multivaries ", \
        "d'apres un nombre de composantes : "
        
    m_estim_nbcomp = v.mixture_estimation(3, 100, [True, True])
    
    m_estim_nbcomp.plot(variable = 1, Title="Estimated mixture")




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
