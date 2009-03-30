""" Convolution """


import interface
import _stat_tool

from _stat_tool import _Convolution
from _stat_tool import _ConvolutionData

__all__ = ['Convolution',
           '_Convolution',
           '_ConvolutionData',]


def Convolution(*args):
    """Construction of an object of type convolution from elementary
    distributions or from an ASCII file.
    
    The distribution of the sum of independent random variables is the
    convolution of the distributions of these elementary random variables.
    
    :Parameters:
      * dist1, dist2, ...(distribution, mixture, convolution, compound) -
        elementary distributions,
      * file_name (string).

    :Returns:    
        If the construction succeeds, the returned object is of type 
        convolution, otherwise no object is returned.

    :Examples:
    .. doctest::
        :options: +SKIP
    
        >>> Convolution(dist1, dist2, ...)
        >>> Convolution(file_name)
      
    .. seealso::    
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`,
        :func:`~openalea.stat_tool.simulate.Simulate`.
    """

    if(len(args)==0):
        raise TypeError()

    # filename
    if(len(args)==1):
        return _stat_tool._Convolution(args[0])

    # build list of distributions
    else:
        return _stat_tool._Convolution(list(args))



# Extend _Convolution
interface.extend_class( _stat_tool._Convolution, interface.StatInterface)


# Extend _ConvolutionData
interface.extend_class( _stat_tool._ConvolutionData, interface.StatInterface)

