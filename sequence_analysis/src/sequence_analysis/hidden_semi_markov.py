#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: hidden_semi_markov.py summary

    A module dedicated to HiddenSemiMarkov objects

    :Code: mature
    :Documentation: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id: hidden_semi_markov.py 9465 2010-08-30 08:46:21Z cokelaer $
    :Usage: >>> from openalea.sequence_analysis import *

.. testsetup::

    from openalea.sequence_analysis import *
"""
__version__ = "$Id: hidden_semi_markov.py 9465 2010-08-30 08:46:21Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _HiddenSemiMarkov

from openalea.stat_tool import error
from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH
from openalea.sequence_analysis._sequence_analysis import OCCUPANCY_THRESHOLD

__all__ = ['HiddenSemiMarkov',
           '_HiddenSemiMarkov']


# Extend dynamically class
interface.extend_class( _HiddenSemiMarkov, interface.StatInterface)

# Add methods to _Vectors


def HiddenSemiMarkov(*args, **kargs):
    """HiddenSemiMarkov

    Construction of an object of type hidden_semi-markov from an ASCII file.

    :Usage:

    .. doctest::
        :options: +SKIP

        >>> HiddenSemiMarkov(file_name, Length=40, Counting=False)

    :Arguments:

    file_name (string)

    :Optional Arguments:

    * Length (int): length of sequences for the computation of the intensity
      and counting characteristic distributions (default value: 20),
    * Counting (bool): computation of counting characteristic distributions
      (default value: True).

    :Returned Object:

    If the construction succeeds, an object of type hidden_semi-markov is
    returned, otherwise no object is returned.

    :Background:

    A hidden semi-Markov chain is constructed from an underlying semi-Markov
    chain (first-order Markov chain representing transition between distinct
    states and state occupancy distributions associated to the nonabsorbing
    states) and nonparametric observation (or state-dependent) distributions.
    The state occupancy distributions are defined as object of type distribution
    with the additional constraint that the minimum time spent in a given
    state is 1 (inf_bound <= 1).

    .. seealso::

        :func:`~openalea.sequence_analysis.compare.Compare` (Markovian models for sequences),
        :func:`~openalea.sequence_analysis.estimate.Estimate` (Markovian models),
        :func:`~openalea.sequence_analysis.data_transform.ComputeStateSequences`,
        :func:`~openalea.sequence_analysis.simulate.Simulate` (Markovian models).
    """
    Length = kargs.get("Length", DEFAULT_LENGTH)
    CountingFlag = kargs.get("Counting", True)
    OldFormat = error.ParseKargs(kargs, "Format", "Current",
                                 {"Current":False, "Old":True})
    error.CheckArgumentsLength(args, 1, 1)
    error.CheckType([args[0], Length, CountingFlag, OldFormat],
                    [str, int, bool, bool])

    filename = args[0]
    if os.path.isfile(filename):

        hsm = _HiddenSemiMarkov(filename, Length, CountingFlag,
                                   OCCUPANCY_THRESHOLD, OldFormat)
    else:
        raise IOError("bad file name %s" % filename)


    return hsm





