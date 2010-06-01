#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: hidden_variable_order_markov.py summary

    A module dedicated to HiddenVariableOrderMarkov objects

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id$
"""
__version__ = "$Id$"

import os

import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import \
    _HiddenVariableOrderMarkov

from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH
from openalea.sequence_analysis._sequence_analysis import OCCUPANCY_THRESHOLD

from openalea.stat_tool import error

__all__ = ['HiddenVariableOrderMarkov',
           '_HiddenVariableOrderMarkov']


# Extend dynamically class
interface.extend_class( _HiddenVariableOrderMarkov, interface.StatInterface)

# Add methods to _Vectors


def HiddenVariableOrderMarkov(filename=None, Length=DEFAULT_LENGTH,
                              CumulThreshold=OCCUPANCY_THRESHOLD):
    """HiddenVariableOrderMarkov
   
    .. todo:: documentation
    """ 
        
    error.CheckType([filename, Length, CumulThreshold], [str, int, float])

    if not os.path.isfile(filename):
        raise IOError("Invalid filename")
    else: 
        return _HiddenVariableOrderMarkov(filename, Length,  CumulThreshold)
        





