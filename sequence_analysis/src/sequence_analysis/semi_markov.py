#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Semi markov module

.. module:: semi_markov
    :synopsis: a module dedicated to SemiMarkov objects

.. topic:: Summary

    A module dedicated to SemiMarkov objects

    :Code: mature
    :Documentation: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id: semi_markov.py 9478 2010-08-31 16:33:16Z cokelaer $
    :Usage: >>> from openalea.sequence_analysis.semi_markov import * 


.. testsetup::
    
    from openalea.sequence_analysis.semi_markov import * 
"""
__version__ = "$Id: semi_markov.py 9478 2010-08-31 16:33:16Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _SemiMarkov
from openalea.sequence_analysis._sequence_analysis import _SemiMarkovData
from openalea.stat_tool import error

from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH

__all__ = ['SemiMarkov',
           '_SemiMarkov',
           '_SemiMarkovData']


# Extend dynamically class
interface.extend_class( _SemiMarkov, interface.StatInterface)
interface.extend_class( _SemiMarkovData, interface.StatInterface)
# Add methods to _Vectors


def SemiMarkov(filename=None, length=DEFAULT_LENGTH, counting=True,
               cumul_threshold=0):
    """SemiMarkov constructor

    Construction of a semi-Markov chain from an ASCII file.

    :Usage:

    ::
    
        SemiMarkov(filename, length=40, counting=True, cumul_threshold=40) 

     .. todo:: make the parameter input names consistent over all modules e.g, length and Length
        should be only denoted either length or Length exclusively. For backward 
        compatibility, Length should be used ?

    :Arguments:

    * filename (string).

    :Optional Arguments:

    * Length (int): length of sequences for the computation of the intensity and
      counting characteristic distributions (default value: 20),
    * Counting (bool): computation of counting characteristic distributions default value: True).

    :Returned Object:

    If the construction succeeds, an object of type semi-markov is returned,
    otherwise no object is returned.

    :Background:

    A semi-Markov chain is constructed from a first-order Markov chain
    representing transition between distinct states and state occupancy
    distributions associated to the non-absorbing states. The state occupancy
    distributions are defined as object of type distribution with the
    additional constraint that the minimum time spent in a given state
    is at least 1 (inf_bound >= 1).

    .. seealso::

        :class:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.sequence_analysis.compare.Compare`,
        :func:`~openalea.sequence_analysis.simulate.Simulate`.

    """
    # todo:: cumul_threshold has no default value in cpp code semi_markov.cpp
    # does zero is a good default value ?
    error.CheckType([filename, length, counting, cumul_threshold],
                    [str, int, bool, [int, float]])

    if not os.path.isfile(filename):
        raise IOError("Invalid filename")
    else:
        return _SemiMarkov(filename, length, counting, cumul_threshold)





