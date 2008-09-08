__doc__ = """ Convolution """
__docformat__ = "restructuredtext"

import interface
import _stat_tool

from _stat_tool import _Convolution
from _stat_tool import _ConvolutionData

__all__ = ['Convolution',
           '_Convolution',
           '_ConvolutionData',
           ]



def Convolution(*args):
    """    
    Construction of an object of type convolution from elementary distributions 
    or from an ASCII file.

    Usage
    -----
      * ``Convolution(dist1, dist2,...)``
      * ``Convolution(file_name)``

    Parameters
    ----------
      * dist1, dist2, ...(distribution, mixture, convolution, compound): elementary distributions,
      * file_name (string). 

    Return
    ------
    If the construction succeeds, the returned object is of type convolution, 
    otherwise no object is returned. 

    Background
    ----------
    The distribution of the sum of independent random variables is the convolution 
    of the distributions of these elementary random variables. 

    See Also
    --------
    ``Save``, ``Estimate``, ``Simulate``.
    """

    if(len(args)==0) : 
        raise TypeError()

    # filename
    if(len(args)==1) :
        return _stat_tool._Convolution(args[0])

    # build list of distributions
    else:
        return _stat_tool._Convolution(list(args))
    


# Extend _Convolution
interface.extend_class( _stat_tool._Convolution, interface.StatInterface)


# Extend _ConvolutionData
interface.extend_class( _stat_tool._ConvolutionData, interface.StatInterface)



########################## Test Convolution ########################################
from openalea.stat_tool import get_test_file

class Test:
    def test_emty(self):
        try:
            m = Convolution()
            assert False

        except Exception:
            assert True


    def test_file(self):
        c = Convolution(get_test_file("convolution1.conv"))
        assert c


    def test_build_convolution(self):
        from distribution import Binomial, NegativeBinomial

        d1 = Binomial(0, 10, 0.5)
        d2 = NegativeBinomial(0, 1, 0.1)

        m = Convolution(d1, d2)
        assert m
        return m


    def test_plot(self):

        m = self.test_build_convolution()
        m.plot()

        assert str(m)
        m.display()


    def test_simulation(self):

        m = self.test_build_convolution()
        s = m.simulate(1000)

        assert len(s) == 1000
        assert s.nb_histogram() == 2
        assert str(s)


    def test_extract(self):
        from data_transform import ExtractDistribution
        from distribution import Binomial, NegativeBinomial

        m = self.test_build_convolution()
        assert m.nb_distribution() == 2

        assert m.extract_convolution() == ExtractDistribution(m, "Convolution")

        assert m.extract_elementary(1) == Binomial(0, 10, 0.5)
        assert m.extract_elementary(2) == NegativeBinomial(0, 1, 0.1)

        assert ExtractDistribution(m, "Elementary", 1) == Binomial(0, 10, 0.5)
        assert ExtractDistribution(m, "Elementary", 2) == NegativeBinomial(0, 1, 0.1)


    def test_extract_data(self):
        from distribution import Binomial

        m = self.test_build_convolution()
        s = m.simulate(1000)

        m = s.estimate_convolution(Binomial(0, 10, 0.5), Estimator="Parametric")

        d = m.extract_data()    
        assert d


