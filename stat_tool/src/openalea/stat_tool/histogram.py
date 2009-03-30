__doc__ = """Histogram functions and classes"""
__docformat__ = "restructuredtext"

import sys
import os
sys.path.append(os.path.abspath("."))

import _stat_tool
import interface

from _stat_tool import _DistributionData

# Extend _DistributionData class dynamically
interface.extend_class(_stat_tool._DistributionData, interface.StatInterface)

__all__ = ["_DistributionData",
           "Histogram",
           ]


def Histogram(arg):
    """Construction of a frequency distribution from an object of type
    list(int) or from an ASCII file.

    In the file syntax, the frequencies *fi* for each possible value
    *i* are given in two columns. In the case of an argument of type
    (list(int)), it is simply assumed that each array element represents
     one data item.

    :Parameters:
      * `list` (list(int)) -
      * `filename` (string) -

    :Returns:
       If the construction succeeds, an object of type `_DistributionData`
       is returned.

    :Examples:
    .. doctest::
        :options: +SKIP

        >>> Histogram(list)
        >>> Histogram(filename)

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

    return _DistributionData(arg)



        


