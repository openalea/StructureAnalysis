"""Compound module

:Status: 
 * constructors done
 * test done
 * error management done
 * documentation to be checked

"""
__revision__ = "$Id$"

import interface
import error

from _stat_tool import _Compound
from _stat_tool import _CompoundData
from _stat_tool import _ParametricModel

__all__ = ['Compound',
            '_Compound',
            '_CompoundData',
            ]


def Compound(*args, **kargs):
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

    .. doctest::
        :options: +SKIP

        >>> Compound(sum_dist, dist)
        >>> Compound(sum_dist, dist, Threshold=0.999)
        >>> Compound(filename)

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`,
        :func:`~openalea.stat_tool.simulate.Simulate`
    """
    error.CheckArgumentsLength(args, 1, 2)
    error.CheckOptionalArgumentsLength(kargs, 0, 1)
    
    Threshold = kargs.get("Threshold", None)

    # filename
    if (len(args)==1):
        error.CheckType(args[0], str, arg_id=1)
        result =  _Compound(args[0])
        
    # build list of distributions
    if (len(args)==2):
        error.CheckType(args[0], _ParametricModel, arg_id=1)
        error.CheckType(args[1], _ParametricModel, arg_id=2)
        if Threshold:
            result =  _Compound(args[0], args[1]) 
        else:
            result =  _Compound(args[0], args[1])
            
    if result is not None:
        return result
    else:
        error.StatToolError('Unknown error in Compound')
    
    
    

# Extend _Compound
interface.extend_class(_Compound, interface.StatInterface)


# Extend _CompoundData
interface.extend_class(_CompoundData, interface.StatInterface)




