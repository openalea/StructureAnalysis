#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""Plot functions

.. topic:: plot.py summary

    A module that provides plotting functions for Gnuplot
    of Matplorlib.

    :Code status: mature
    :Documentation status: to be completed
    :Authors: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>,
        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>

    :Revision: $Id$

.. inheritance-diagram:: openalea.stat_tool.plot

"""
__version__ = "$Id$"

__all__ = []

from functools import wraps

from __stat_tool.stat_tool import SinglePlot, MultiPlot, MultiPlotSet
from std import Iterator

def __iter__(self):
    return Iterator(self.begin(), self.end())

SinglePlot.__iter__ = __iter__
del __iter__

#MultiPlot.__len__ = MultiPlot.size
#del MultiPlot.size

def wrapper(f):
    @wraps(f)
    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError()
        return f(self, index)
    return __getitem__

MultiPlot.__getitem__ = wrapper(MultiPlot.__getitem__)
del wrapper #, MultiPlot.__getitem__

#MultiPlotSet.__len__ = MultiPlotSet.size
#del MultiPlotSet.size

def wrapper(f):
    @wraps(f)
    def __getitem__(self, index):
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError()
        return f(self, index)
    return __getitem__

MultiPlotSet.__getitem__ = wrapper(MultiPlotSet.__getitem__)
del wrapper, MultiPlotSet.__getitem__