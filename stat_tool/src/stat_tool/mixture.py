#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Mixture object

.. topic:: mixture.py summary

    A module dedicated to Mixture objects

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: mixture.py 15183 2013-11-06 10:35:50Z jbdurand $
    
.. testsetup:: *

    from openalea.stat_tool.vectors import *
"""
__version__ = "$Id: mixture.py 15183 2013-11-06 10:35:50Z jbdurand $"


import interface
import error

from openalea.stat_tool._stat_tool import _DiscreteMixture
from openalea.stat_tool._stat_tool import _DiscreteMixtureData
from openalea.stat_tool._stat_tool import _DiscreteParametricModel
from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _Convolution

__all__ = ['_DiscreteMixture',
           '_DiscreteMixtureData',
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
      * `dist1`, `dist2`, ... (`_DiscreteParametricModel`, `_DiscreteMixture`, `_Convolution`,
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

    types = [_DiscreteParametricModel, _DiscreteMixture, _Compound, _Convolution]

    # filename 
    if (len(args) == 1):
        error.CheckType([args[0]], [str], arg_id=[1])
        result = _DiscreteMixture(args[0])

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

        result = _DiscreteMixture(weights, dists)

    return result

# Extend _DiscreteMixture
interface.extend_class(_DiscreteMixture, interface.StatInterface)

# Extend _DiscreteMixtureData
interface.extend_class(_DiscreteMixtureData, interface.StatInterface)

_DiscreteMixture.__doc__ = Mixture.__doc__
