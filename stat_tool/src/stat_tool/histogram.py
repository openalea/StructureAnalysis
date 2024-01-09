#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Histogram module

.. topic:: histogram.py summary

    A module dedicated to Histogram and DiscreteDistribution

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: histogram.py 9453 2010-08-28 08:10:59Z cokelaer $
"""
__version__ = "$Id: histogram.py 9453 2010-08-28 08:10:59Z cokelaer $"

#import sys
#import os
#sys.path.append(os.path.abspath("."))

import interface

from openalea.stat_tool._stat_tool import _DiscreteDistributionData

# Extend _DistributionData class dynamically
interface.extend_class(_DiscreteDistributionData, interface.StatInterface)

__all__ = ["_DiscreteDistributionData",
           "Histogram",
           ]


def Histogram(*args):
    """Construction of a frequency distribution from an object of type
    list(int) or from an ASCII file.

    In the file syntax, the frequencies *fi* for each possible value
    *i* are given in two columns. In the case of an argument of type
    (list(int)), it is simply assumed that each array element represents
    one data item.

    :param integer list: a list of integer values
    :param string filename: a valid filename in the proper format (see syntax part)

    :Returns:
       If the construction succeeds, an object of type `_DiscreteDistributionData`
       is returned.

    :Usage:

    .. doctest::
        :options: +SKIP

        >>> Histogram(list)
        >>> Histogram(filename)

    :Examples:

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool import *
        from numpy.random import randint
        h = Histogram(randint(10, size=10000).tolist())
        h.plot()

    .. note:: works for integer values only. 

    .. seealso::
        :func:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.stat_tool.cluster.Cluster`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.stat_tool.data_transform.Shift`,
        :func:`~openalea.stat_tool.cluster.Transcode`,
        :func:`~openalea.stat_tool.data_transform.ValueSelect`,
        :func:`~openalea.stat_tool.comparison.Compare`,
        :func:`~openalea.stat_tool.estimate.Estimate`

    """

    ret = None
    # Histogram(filename)
    if len(args)==1 and isinstance(args[0], str):
        try:
            ret = _DiscreteDistributionData(args[0])
        except:
            raise IOError("wrong filename ? %s" % args[0])
    # Histogram([1,2,3])
    elif len(args)==1 and isinstance(args[0], list):
        ret = _DiscreteDistributionData(args[0])
    # Histogram(1,2,3)
    else:
        ret = _DiscreteDistributionData(list(args))

    return ret
