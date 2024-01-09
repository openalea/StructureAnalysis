#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""time_events module

.. module:: time_events

.. topic::  summary

    A module dedicated to TimeEvents objects

    :Code: mature
    :Documentation: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: time_events.py 15827 2014-02-19 16:17:45Z jbdurand $
    
"""
__version__ = "$Id: time_events.py 15827 2014-02-19 16:17:45Z jbdurand $"

import os
import openalea.stat_tool.interface as interface
from openalea.stat_tool._stat_tool import _DiscreteMixtureData
from openalea.stat_tool._stat_tool import _CompoundData
from openalea.stat_tool._stat_tool import _ConvolutionData
from openalea.stat_tool._stat_tool import _DiscreteDistributionData

from openalea.sequence_analysis._sequence_analysis import _TimeEvents
from openalea.sequence_analysis._sequence_analysis import _RenewalData
from openalea.sequence_analysis._sequence_analysis import _Sequences

#import _sequence_analysis
from openalea.stat_tool import error

__all__ = ['TimeEvents',
           '_TimeEvents',
           'NbEventSelect']



# Extend dynamically class
interface.extend_class( _TimeEvents, interface.StatInterface)

# Add methods to _Vectors


def TimeEvents(*args, **kargs):
    """TimeEvents

    Construction of data of type {time interval between two observation dates,
    number of events occurring between these two observation dates} from time
    sequences, from an object of type HISTOGRAM or from an ASCII file.

    :Usage:

    .. doctest::
        :options: +SKIP
        
        >>> TimeEvents(seq1, begin_date, end_date, PreviousDate=3, NextDate=8)
        >>> TimeEvents(seqn, variable, begin_date, end_date, PreviousDate=3,\
                NextDate=8)
        >>> TimeEvents(histo, time)
        >>> TimeEvents(file_name)
        >>> h = Histogram([1,1,1,2,2,2])
        >>> t = TimeEvents(h, 2)

    :Arguments:

    * seq1 (sequences): univariate time sequences (with an explicit index
      parameter of type TIME),
    * seqn (sequences): multivariate time sequences (with an explicit index
      parameter of type TIME),
    * variable (int): variable index,
    * begin_date (int): initial observation date,
    * end_date (int): final observation date,
    * histo (histogram, mixture_data, convolution_data, compound_data): number of
      events frequency distribution,
    * time (int): time interval between two observation dates (length of the
      observation period),
    * file_name (string).

    :Optional Arguments:

    * PreviousDate (int): date preceding the initial observation date to check
      the increasing character of the number of events. This optional argument
      can only be used if the first mandatory argument is of type sequences.
    * NextDate (int): date following the final observation date to check the
      increasing character of the number of events. This optional argument can
      only be used if the first mandatory argument is of type sequences.

    :Returned Object:

    If the construction succeeds, an object of type time_events is returned,
    otherwise no object is returned.

    .. seealso::

        :func:`Save`,
        :func:`~openalea.stat_tool.data_transform.ExtractHistogram`,
        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.sequence_analysis.time_events.NbEventSelect`,
        :func:`~openalea.sequence_analysis.data_transform.TimeScaling`,
        :func:`~openalea.sequence_analysis.data_transform.TimeSelect`.

    .. todo:: fix the build_time_events method to allows constructor with
       histogram issue: this method is in stat_tool and returns a time events
       so stat_tool requires to know sequence_analysis...
    """

    PreviousDate = kargs.get("PreviousDate", -1)
    NextDate = kargs.get("NextDate", -1)

    if len(args) == 1 and isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            time_events =  _TimeEvents(filename)
        else:
            raise IOError("bad file name")
    elif isinstance(args[0], _Sequences):
        seq = args[0]
        nb_variable = seq.nb_variable
        if nb_variable != 1:
            variable = args[0]
            begin_date = args[1]
            end_date = args[2]
        else:
            variable = 1
            begin_date = args[1]
            end_date = args[2]
        error.CheckType([variable, begin_date, end_date], [int, int, int])

        time_events = seq.extract_time_events(variable, begin_date, end_date,
                                     PreviousDate, NextDate)

    else:
        # should work with Histogram, Mixture_data, Conv_data, comp_data
        error.CheckArgumentsLength(args, 2, 2)
        error.CheckType([args[0], args[1]], \
                        [[_DiscreteDistributionData, _DiscreteMixtureData,\
                          _ConvolutionData, _CompoundData], int])
        distribution = args[0]
        time = args[1]
        time_events = _TimeEvents(distribution, time)


    return time_events


def NbEventSelect(obj, imin, imax):
    """NbEventSelect

    Selection of data item of type {time interval between two observation dates, 
    number of events occurring between these two observation dates} according to
    a number of events criterion.

    :Usage:

    .. doctest:: 
        :options: +SKIP
        
        >>> NbEventSelect(timev, min_nb_event, max_nb_event)


    :param TimeEvents,RenewalData time_v:
    :param int min_nb_event: minimum number of events,
    :param int max_nb_event: maximum number of events.

    :Returned Object:

    If 0 <= min_nb_event < max_nb_event and if the range of number of events defined by min_nb_event and max_nb_event enables to select data items of type {time interval between two observation dates, number of events}, an object of type time_events is returned, otherwise no object is returned.

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.Merge`,
        :func:`~openalea.sequence_analysis.data_transform.TimeScaling`,
        :func:`~openalea.sequence_analysis.data_transform.TimeSelect`.
    """
    error.CheckType([obj, imin, imax],
                    [[_TimeEvents, _RenewalData], int, int])

    return obj.nb_event_select(imin, imax)
