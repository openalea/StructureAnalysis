""" Mixture class

:Authors: Thomas Cokelaer <Thomas.Cokelaer@inria.fr> 
 
todo: MixtureData ?
"""
__version__ = "$Id$"


import interface
import error

from _stat_tool import _Mixture
from _stat_tool import _MixtureData
from _stat_tool import _ParametricModel
from _stat_tool import _Compound
from _stat_tool import _Convolution

__all__ = ['_Mixture',
           '_MixtureData',
           'Mixture', ]


def Mixture(*args):
    """Construction of a mixture of distributions from elementary distributions
    and associated weights or from an ASCII file.

    A mixture is a parametric model of classification where each elementary
    distribution or component represents a class with its associated weight.

    :Parameters:
      * `weight1`, `weight2`, ... (float) - weights of each component.
         These weights should sum to one (they constitute a discrete
         distribution).
      * `dist1`, `dist2`, ... (`_ParametricModel`, `_Mixture`, `_Convolution`,
        `_Compound`) elementary distributions (or components).
      * `filename` (string) -

    :Returns:
        If the construction succeeds, an object of type mixture is returned,
        otherwise no object is returned.

    :Examples:

    .. doctest::
        :options: +SKIP

        >>> Mixture(weight1, dist1, weight2, dist2,...)
        >>> Mixture(filename)

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`,
        :func:`~openalea.stat_tool.simulate.Simulate`.

    """
    error.CheckArgumentsLength(args, 1)

    types = [_ParametricModel, _Mixture, _Compound, _Convolution]

    # filename 
    if (len(args) == 1):
        error.CheckType([args[0]], [str], arg_id=[1])
        result = _Mixture(args[0])

    # build list of weights and distributions
    else:
        nb_param = len(args)
        if ((nb_param % 2) != 0):
            raise TypeError("Number of parameters must be pair")

        # Extract weights ands distributions
        weights = []
        dists = []
        for i in xrange(nb_param / 2):
            weights.append(args[i * 2])
            error.CheckType([args[i * 2 + 1]], [types], arg_id=[i * 2 + 1])
            error.CheckType([args[i * 2]], [float], arg_id=[i * 2])
            #dists.append(_Distribution(args[i * 2 + 1]))
            dists.append((args[i * 2 + 1]))

        result = _Mixture(weights, dists)

    return result

# Extend _Mixture
interface.extend_class(_Mixture, interface.StatInterface)

# Extend _MixtureData
interface.extend_class(_MixtureData, interface.StatInterface)


