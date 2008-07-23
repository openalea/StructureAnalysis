__docformat__ = "restructuredtext"
__doc__ = """ Mixture """


import interface
import _stat_tool

from _stat_tool import _Mixture
from _stat_tool import _MixtureData

__all__ = ['_Mixture', 
           '_MixtureData',
           'Mixture',
           ]




def Mixture(*args):
    """
    Construction of a mixture of distributions from elementary distributions and 
    associated weights or from an ASCII file.

    Usage
    -----
      * ``Mixture(weight1, dist1, weight2, dist2,...)``
      * ``Mixture(filename)``

    Parameters
    ----------
      * weight1, weight2, ... (float): weights of each component. \
     These weights should sum to one (they constitute a discrete distribution).
      * dist1, dist2, ... (`_ParametricModel`, `_Mixture`, `_Convolution`, `_Compound`)\ 
      elementary distributions (or components).
      * filename (string). 

    Return
    ------
    If the construction succeeds, an object of type mixture is returned, otherwise no object 
    is returned. 

    Background
    ----------
    A mixture is a parametric model of classification where each elementary distribution or 
    component represents a class with its associated weight. 

    See Also
    --------
    `Save`, `Estimate`, `Simulate`.

    """

    if(len(args)==0) : 
        raise TypeError()

    # filename
    if(len(args)==1) :
        return _stat_tool._Mixture(args[0])

    # build list of weights and distributions
    else:
        nb_param = len(args)

        if((nb_param % 2) != 0) :
            raise TypeError("Number of parameters must be pair")

        weights = []
        dists = []
        
        for i in xrange(nb_param / 2):
            weights.append(args[i * 2])
            dists.append(args[i * 2 + 1])

        return _stat_tool._Mixture(weights, dists)
    


# Extend _Mixture
interface.extend_class( _stat_tool._Mixture, interface.StatInterface)


# Extend _MixtureData
interface.extend_class( _stat_tool._MixtureData, interface.StatInterface)



########################## Test Mixture ########################################

class Test:
    def test_emty(self):

        try:
            m = Mixture()
            assert False

        except Exception:
            assert True


    def test_file(self):

        m = Mixture("../../../test/mixture1.mixt")
        assert m


    def test_build_mixture(self):

        from distribution import Uniform


        d1 = Uniform(0, 10)
        d2 = Uniform(10, 20)
        d3 = Uniform(20, 30)

        m = Mixture(0.1, d1, 0.2, d2, 0.7, d3)
        assert m
        return m


    def test_plot(self):

        m = self.test_build_mixture()
        #    m.plot()

        assert str(m)
        m.display()


    def test_simulation(self):

        m = self.test_build_mixture()
        s = m.simulate(10)

        assert s.nb_component() == 3
        assert str(s)


    def test_extract(self):

        from data_transform import ExtractDistribution
        from distribution import Uniform


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

        from histogram import Histogram 

        h = Histogram("../../../test/meri2.his")
        m = h.estimate_mixture(["B", "NB"])

        d = m.extract_data()
        assert d


