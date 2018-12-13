#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Histogram module

.. topic:: histogram.py summary

    A module dedicated to Histogram and DiscreteDistribution

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

"""

from . import error
from .enum import *

from stat_tool.__stat_tool.stat_tool import *
import stat_tool.__stat_tool.stat_tool as cst

_DiscreteDistributionData = cst.DiscreteDistributionData

__all__ = ["_DiscreteDistributionData",
           "FrequencyDistribution",
           "Histogram",
           ]


def FrequencyDistribution(*args):
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

        >>> FrequencyDistribution(list)
        >>> FrequencyDistribution(filename)

    :Examples:

    .. plot::
        :width: 50%
        :include-source:

        from openalea.stat_tool import *
        from numpy.random import randint
        h = FrequencyDistribution(randint(10, size=10000).tolist())
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
            #ret = _DiscreteDistributionData(args[0])
            ret =  _DiscreteDistributionData.ascii_read(args[0])

        except:
            raise IOError("wrong filename ? %s" % args[0])
    # Histogram([1,2,3])
    elif len(args)==1 and isinstance(args[0], list):
        n = len(args[0])
        ret = _DiscreteDistributionData(n, args[0])
    # Histogram(1,2,3)
    else:
        l = list(args)
        n = len(l)
        ret = _DiscreteDistributionData(n, l)

    return ret


Histogram = FrequencyDistribution
