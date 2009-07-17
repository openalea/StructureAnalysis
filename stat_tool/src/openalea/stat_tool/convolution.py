""" Convolution module 

:Status:
  * constructors done
  * test done
  * error management done
  * documentation to be checked
  
"""
__revision__ = "$Id$"


import interface
import error

#import _stat_tool

from _stat_tool import _Convolution
from _stat_tool import _ConvolutionData
from _stat_tool import _ParametricModel

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
    error.CheckArgumentsLength(args, 1)

    # filename
    if(len(args)==1):
        error.CheckType(args[0], str, arg_id=1)
        result =  _Convolution(args[0])
    # build list of distributions
    else:
        #check that all arguments are _parametric_model
        for arg, i in zip(args, range(0, len(args))):
            print 'check %s' %i
            error.CheckType(arg, _ParametricModel, arg_id=i+1)
        
        result = _Convolution(list(args))
        
    if result is not None:
        return result
    else:
        error.StatToolError('Unknown error in Convolution')



# Extend _Convolution
interface.extend_class(_Convolution, interface.StatInterface)


# Extend _ConvolutionData
interface.extend_class(_ConvolutionData, interface.StatInterface)

