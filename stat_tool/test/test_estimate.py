"""estimate tests"""
__revision__ = "$Id: $"

from openalea.stat_tool import distribution, data_transform, histogram
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.data_transform import Shift
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.estimate import Estimate
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool import  *


class Test:

    def test_nonparametric(self):

        h = Histogram(("meri1.his"))
        e =  h.estimate_nonparametric()
        assert e
        # assert abs(e.likelihood() - VAL) < epsilon

    def test_nb(self):
        """NegativeBinomial"""
        h = Histogram(("peup2.his"))
        assert h.estimate_parametric('NB')


    def test_binomial(self):
        """BINOMIAL Distribution"""
        h = Histogram(("meri5.his"))
        assert h.estimate_parametric('B')

    def test_poisson(self):
        """Poisson distribution
        >>> p = distribution.Poisson(0, 10)
        >>> h = p.simulate(1000)        
        """
        p = distribution.Poisson(0, 10)
        h = p.simulate(1000)

        assert h.estimate_parametric('P')


    def test_mixture_1(self):

        h = Histogram(("peup2.his"))
        m1 =  h.estimate_mixture(["B", "NB", "NB", "NB"], NbComponent="Estimated")
        assert m1


    def test_mixture_2(self):

        h = Histogram(("peup2.his"))
        m2 = h.estimate_mixture([Binomial(0, 10, 0.5), "NB"])
        assert m2


    def test_convolution(self):

        elementary = Histogram(("nothofagus_antarctica_bud_2.his"))
        total = Histogram(("nothofagus_antarctica_shoot_2.his"))

        convol1 = Estimate(Shift(total, 1), "CONVOLUTION", Estimate(elementary, "NP"), 
                           NbIteration=100, Estimator="PenalizedLikelihood", Weight=0.5)

        convol2 = total.shift(1).estimate_convolution(elementary.estimate_nonparametric(), 
                                                      NbIteration=100, 
                                                      Estimator="PenalizedLikelihood", 
                                                      Weight=0.5)

        assert convol1 and convol2
        assert convol1 == convol2



    def test_compound(self):
        
        cdist1 = Compound("compound1.cd")
        chisto1 = Simulate(cdist1, 200)
        cdist2 = Estimate(chisto1, "COMPOUND",
                      ExtractDistribution(cdist1, "Elementary"),
                      ExtractDistribution(cdist1, "Sum"),
                      MinInfBound=0)
    
        cdist3 = chisto1.estimate_compound(
                                  ExtractDistribution(cdist1, "Elementary"), 
                                  ExtractDistribution(cdist1, "Sum"))
        
        assert cdist2==cdist3
        


 


