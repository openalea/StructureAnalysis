#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: nonhomogeneous_markov.py summary

    A module dedicated to NonhomogeneousMarkov objects

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id$
"""
__revision__ = "$Id$"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _NonHomogeneousMarkov
from openalea.sequence_analysis._sequence_analysis import _NonHomogeneousMarkovData
from openalea.stat_tool import error
from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH


__all__ = ['NonhomogeneousMarkov',
           '_NonHomogeneousMarkov',
           '_NonHomogeneousMarkovData',
           ]


# Extend dynamically class
interface.extend_class( _NonHomogeneousMarkov, interface.StatInterface)
interface.extend_class( _NonHomogeneousMarkovData, interface.StatInterface)


def NonhomogeneousMarkov(filename, length=DEFAULT_LENGTH):
    """NonhomogeneousMarkov constructor

    :param filename:
    :param length: optional argument (default is 20)

    :Usage:
    
    .. doctest:: 
        :options: +SKIP
        
        >>> nm = NonhomogeneousMarkov("filename.dat")
        >>> nm = NonhomogeneousMarkov("filename.dat", 10)
    """
    error.CheckType([filename, length], [str, int])

    if os.path.isfile(filename):
        return _NonHomogeneousMarkov(filename, length)
    else:
        raise IOError("bad file name")





