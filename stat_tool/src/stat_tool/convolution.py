#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Convolution module

.. topic:: Summary

    A module dedicated to Convolution

    :Code status: mature
    :Documentation status: mature
    :Author: Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
        Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: convolution.py 15183 2013-11-06 10:35:50Z jbdurand $

.. testsetup:: *

    from openalea.stat_tool import *
"""
__version__ = "$Id: convolution.py 15183 2013-11-06 10:35:50Z jbdurand $"


import interface
import error


from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _DiscreteMixture
from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _ConvolutionData
from openalea.stat_tool._stat_tool import _DiscreteParametricModel

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

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool import *
        sum_dist = Binomial(0,10,0.5)
        dist = Binomial(0,15,0.2)
        c = Convolution(sum_dist, dist)
        c.plot()


    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.estimate.Estimate`,
        :func:`~openalea.stat_tool.simulate.Simulate`.
    """
    error.CheckArgumentsLength(args, 1)

    possible_types = [_DiscreteParametricModel, _DiscreteMixture,
                      _Compound, _Convolution]

    # filename
    if(len(args)==1):
        error.CheckType([args[0]], [str], arg_id=[1])
        result =  _Convolution(args[0])
    # build from list of distributions
    else:
        arguments = []
        #check that all arguments are correct
        for arg, i in zip(args, range(0, len(args))):
            error.CheckType([arg], [possible_types], variable_pos=[i+1])
            arguments.append(arg)
        result = _Convolution(arguments)

    return result


# Extend _Convolution
interface.extend_class(_Convolution, interface.StatInterface)


# Extend _ConvolutionData
interface.extend_class(_ConvolutionData, interface.StatInterface)

