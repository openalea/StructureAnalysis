#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""nonhomogeneous_markov module

.. module:: nonhomogeneous_markov
    :synopsis: A module dedicated to NonhomogeneousMarkov objects

.. topic:: nonhomogeneous_markov.py summary

    A module dedicated to NonhomogeneousMarkov objects

    :Code: mature
    :Documentation: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>
    :Revision: $Id$
    :Usage: from openalea.sequence_analysis.nonhomogeneous_markov import *

"""
__revision__ = "$Id$"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _NonHomogeneousMarkov, _NonHomogeneousMarkovData, DEFAULT_LENGTH
from openalea.stat_tool import error


__all__ = ['NonhomogeneousMarkov',
           '_NonHomogeneousMarkov',
           '_NonHomogeneousMarkovData',
           ]


# Extend dynamically class
interface.extend_class( _NonHomogeneousMarkov, interface.StatInterface)
interface.extend_class( _NonHomogeneousMarkovData, interface.StatInterface)


def NonhomogeneousMarkov(filename, length=DEFAULT_LENGTH):
    """NonhomogeneousMarkov constructor

    :param str filename:
    :param int length: 

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





