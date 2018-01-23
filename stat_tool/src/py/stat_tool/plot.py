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

from . import _stat_tool
from .__stat_tool.stat_tool import SinglePlot, MultiPlot, MultiPlotSet

class Point(object):

    def __init__(self, single_plot, index):
        if index < 0:
            index += len(single_plot)
        if not 0 <= index < len(single_plot):
            raise IndexError()
        self._single_plot = single_plot
        self._index = index

    @property
    def first(self):
        return self._single_plot.get_x(self._index)

    @property
    def second(self):
        return self._single_plot.get_y(self._index)

    @property
    def label(self):
        return self._single_plot.get_label(self._index)

def __getitem__(self, index):
    return Point(self, index)

SinglePlot.__getitem__ = __getitem__
del __getitem__

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
del wrapper

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
del wrapper