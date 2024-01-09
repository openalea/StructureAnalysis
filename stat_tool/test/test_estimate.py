"""estimate tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__version__ = "$Id: test_estimate.py 9397 2010-08-10 10:50:07Z cokelaer $"


from openalea.stat_tool import distribution
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.data_transform import Shift, ExtractDistribution
from openalea.stat_tool.histogram import Histogram

from openalea.stat_tool.compound import Compound
from openalea.stat_tool.simulate import Simulate
from openalea.stat_tool import  _stat_tool
from openalea.stat_tool.enums import distribution_identifier_type

from openalea.stat_tool.estimate import Estimate, likelihood_penalty_type

from tools import runTestClass


from openalea.stat_tool import get_shared_data

class Test():

    def __init__(self):
        pass

    def test_nonparametric(self):
        h = Histogram(get_shared_data("meri1.his"))
        e =  h.estimate_nonparametric()
        assert e
        # assert abs(e.likelihood() - VAL) < epsilon

    def test_nb(self):
        """NegativeBinomial"""
        h = Histogram(get_shared_data("peup2.his"))
        assert h.estimate_parametric('NB')

    def test_binomial(self):
        """BINOMIAL Distribution"""
        h = Histogram(get_shared_data("meri5.his"))
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
        distributions = ["B", "NB", "NB", "NB"]
        h = Histogram(get_shared_data( "peup2.his"))
        m1 =  h.estimate_mixture(distributions, NbComponent="Estimated")
        assert m1

        types = []
        for d in distributions:
            temp = distribution_identifier_type[d]
            types.append(temp)

        c = h.mixture_estimation2(types, 0, True, True,
                                likelihood_penalty_type['AIC'])

        assert str(c)==str(m1)

    def test_mixture_2(self):
        h = Histogram(get_shared_data( "peup2.his"))
        m2 = h.estimate_mixture([Binomial(0, 10, 0.5), "NB"])
        assert m2

    def test_convolution(self):

        elementary = Histogram("data/nothofagus_antarctica_bud_2.his")
        total = Histogram("data/nothofagus_antarctica_shoot_2.his")

        convol1 = Estimate(Shift(total, 1), "CONVOLUTION",
                           Estimate(elementary, "NP"),
                           NbIteration=100,
                           Estimator="PenalizedLikelihood",
                           Weight=0.5)

        convol2 = total.shift(1).estimate_convolution(
                                    elementary.estimate_nonparametric(),
                                    NbIteration=100,
                                    Estimator="PenalizedLikelihood",
                                    Weight=0.5)

        assert convol1 and convol2
        assert convol1 == convol2

    def test_compound_two_distribution(self):
        """test to be checked"""
        cdist1 = Compound("data/compound1.cd")
        chisto1 = Simulate(cdist1, 200)

        cdist2 = Estimate(chisto1, "COMPOUND",
                      ExtractDistribution(cdist1, "Sum"),
                      "Sum",
                      InitialDistribution=\
                        ExtractDistribution(cdist1, "Elementary"))

        # If we call the method directly, we need to provide
        # the default values and perform a conversion.
        # Default is LIKELIHOOD -1 -1.0 SECOND_DIFFERENCE ZERO, which
        #  corresponds to 0, -1,-1,1,0
        # In addition because the type is 's', the 2 distributions
        # must be reversed.

        cdist3 = chisto1.compound_estimation1(
                    ExtractDistribution(cdist1, "Elementary"),
                    ExtractDistribution(cdist1, "Sum"),
                    's', _stat_tool.LIKELIHOOD, -1, -1.,
                    _stat_tool.SECOND_DIFFERENCE,
                    _stat_tool.ZERO)
        assert str(cdist2) == str(cdist3)




if __name__ == "__main__":
    runTestClass(Test())

