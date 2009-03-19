
"""Compound module"""
__revision__ = "$Id$"

import interface
import _stat_tool


from _stat_tool import _Compound
from _stat_tool import _CompoundData

__all__ = ['Compound',
            '_Compound',
            '_CompoundData',
            ]


def Compound(*args):
    """
    Construction of a compound of distributions from a sum distribution and an
    elementary distribution or from an ASCII file.
    
    A compound (or stopped-sum) distribution is defined as the distribution 
    of the sum of n independent and identically distributed random variables :math:`X_i`
    where `n` is the value taken by the random variable `N`. The distribution of N is referred
    to as the sum distribution while the distribution of the :math:`X_i` is referred to as
    the elementary distribution.

    :Parameters:
      * `sum_dist` (DISTRIBUTION, MIXTURE, CONVOLUTION, COMPOUND) - sum distribution
      * `dist` (DISTRIBUTION, MIXTURE, CONVOLUTION, COMPOUND) - elementary distribution
      * `filename` (STRING) -

    :Returns:
    
        If the construction succeeds, an object of type `COMPOUND` is returned,
        otherwise no object is returned.
        
    :Examples:
        >>> Compound(sum_dist, dist)
        >>> Compound(filename)

    .. seealso::
        `Save`, `Estimate`, `Simulate`
    """
    if((len(args)==0) or (len(args)>3)) : 
        raise TypeError("Bad number of arguments")

    # filename
    if(len(args)==1) :
        return _stat_tool._Compound(args[0])

    # build list of distributions
    if(len(args)==2) :    
        return _stat_tool._Compound(args[0], args[1])
    
    if(len(args)==3) :    
        return _stat_tool._Compound(args[0], args[1], args[2])
    
    


# Extend _Compound
interface.extend_class(_Compound, interface.StatInterface)


# Extend _CompoundData
interface.extend_class(_CompoundData, interface.StatInterface)



######################## Test Compound ########################################
from openalea.stat_tool import get_test_file


class Test:

    def test_emty(self):
        try:
            m = Compound()
            assert False

        except TypeError:
            assert True

    def test_file(self):
        c = Compound(get_test_file("compound1.cd"))
        assert c

    def test_build_compound(self):
        from distribution import Uniform

        d1 = Uniform(0,10)
        d2 = Uniform(10,20)

        m = Compound(d1, d2)
        assert m
        return m

    def test_plot(self):

        m = self.test_build_compound()
        m.plot()

        assert str(m)
        m.display()

    def test_simulation(self):

        m = self.test_build_compound()
        s = m.simulate(1000)

        assert len(s) == 1000
        assert str(s)

    def test_extract(self):
        from data_transform import ExtractDistribution
        from distribution import Uniform

        m = self.test_build_compound()

        assert m.extract_compound() == ExtractDistribution(m, "Compound")

        assert m.extract_sum() == Uniform(0,10)
        assert m.extract_sum() == ExtractDistribution(m, "Sum")

        assert m.extract_elementary() == Uniform(10,20)
        assert m.extract_elementary() == ExtractDistribution(m, "Elementary")

    def _test_extract_data(self):
        # to be done
        assert False


