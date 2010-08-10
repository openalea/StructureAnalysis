#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Compound object

.. topic:: compound.py summary

    A module dedicated to Compound objects

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id$

"""
__version__ = "$Id$"

import interface
import error

from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _Mixture
from openalea.stat_tool._stat_tool import _CompoundData
from openalea.stat_tool._stat_tool import _DiscreteParametricModel

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
    error.CheckKargs(kargs, possible_kargs = ["Threshold"])

    Threshold = kargs.get("Threshold", None)

    # filename
    if len(args)==1:
        error.CheckType([args[0]], [str])
        result =  _Compound(args[0])

    possible_types = [_DiscreteParametricModel, _Mixture,
                      _Compound, _Convolution]

    # build from two objects and optional threshold
    if len(args)==2:
        error.CheckType([args[0], args[1]],
                        [possible_types, possible_types],
                        variable_pos=[1,2])

        if Threshold:
            result =  _Compound([args[0], args[1]], Threshold)
        else:
            result =  _Compound([args[0], args[1]])

    return result


# Extend _Compound
interface.extend_class(_Compound, interface.StatInterface)


# Extend _CompoundData
interface.extend_class(_CompoundData, interface.StatInterface)




