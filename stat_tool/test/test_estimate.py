"""estimate tests"""
__revision__ = "$Id: $"

from openalea.stat_tool import distribution, data_transform, histogram
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.data_transform import Shift
from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.estimate import Estimate


class Test:

    def test_nonparametric(self):

        h = Histogram(("meri1.his"))
        e =  h.estimate_nonparametric()
        assert e

        # assert abs(e.likelihood() - VAL) < epsilon


    def test_nb(self):

        h = Histogram(("peup2.his"))
        assert h.estimate_parametric('NB')


    def test_binomial(self):

        h = Histogram(("meri5.his"))
        assert h.estimate_parametric('B')


    def test_poisson(self):
        """        
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



    def _test_compound(self):

        assert False





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


