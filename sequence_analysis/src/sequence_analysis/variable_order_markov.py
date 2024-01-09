#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""variable order markov

.. module:: variable_order_markov
    :synopsis: a module dedicated to Variable Order Markov objects

.. topic:: variable_order_markov.py summary

    A module dedicated to Variable Order Markov objects

    :Code: mature
    :Documentation: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id: variable_order_markov.py 9885 2010-11-06 18:19:34Z cokelaer $
    :Usage: >>> from openalea.sequence_analysis import *

"""
__version__ = "$Id: variable_order_markov.py 9885 2010-11-06 18:19:34Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import \
    _VariableOrderMarkov, _VariableOrderMarkovData
from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH

from openalea.stat_tool import error

__all__ = ['VariableOrderMarkov',
           '_VariableOrderMarkov',
           '_VariableOrderMarkovData']


# Extend dynamically class
interface.extend_class( _VariableOrderMarkov, interface.StatInterface)
interface.extend_class( _VariableOrderMarkovData, interface.StatInterface)

# Add methods to _Vectors


def VariableOrderMarkov(*args, **kargs):
    """VariableOrderMarkov

    :Usage:

    .. doctest::
        :options: +SKIP

        >>> VariableOrderMarkov(filename)
    """
    error.CheckArgumentsLength(args, 1, 1)
    error.CheckType([args[0]], [str])
    filename = args[0]
    Length = kargs.get("Length", DEFAULT_LENGTH)

    if os.path.isfile(filename):
        vom =  _VariableOrderMarkov(filename, Length)
    else:
        raise IOError("bad file name %s" % filename)


    return vom




