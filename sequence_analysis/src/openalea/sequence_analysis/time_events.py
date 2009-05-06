"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Time_events

import _sequence_analysis

__all__ = ['TimeEvents',
           '_Time_events']


# Extend dynamically class
interface.extend_class( _Time_events, interface.StatInterface)

# Add methods to _Vectors


def TimeEvents(filename=None):
    """TimeEvents
    
    Construction of data of type {time interval between two observation dates,
     number of events occurring between these two observation dates} from time 
     sequences, from an object of type HISTOGRAM or from an ASCII file.
     
    :Usage:
    
    >>> TimeEvents(seq1, begin_date, end_date, PreviousDate->3, NextDate->8)   
    >>> TimeEvents(seqn, variable, begin_date, end_date, PreviousDate->3, NextDate->8)   
    >>> TimeEvents(histo, time)
    >>> TimeEvents(file_name)
        
    :Arguments:
    
    * seq1 (sequences): univariate time sequences (with an explicit index parameter of type TIME),
    * seqn (sequences): multivariate time sequences (with an explicit index parameter of type TIME),
    * variable (int): variable index,
    * begin_date (int): initial observation date,
    * end_date (int): final observation date,
    * histo (histogram, mixture_data, convolution_data, compound_data): number of 
      events frequency distribution,
    * time (int): time interval between two observation dates (length of the 
      observation period),
    * file_name (string).
    
    :Optional Arguments: 
    
    PreviousDate (int): date preceding the initial observation date to check 
    the increasing character of the number of events. This optional argument 
    can only be used if the first mandatory argument is of type sequences.
    NextDate (int): date following the final observation date to check the 
    increasing character of the number of events. This optional argument can 
    only be used if the first mandatory argument is of type sequences.
    
    :Returned Object:
    
    If the construction succeeds, an object of type time_events is returned, 
    otherwise no object is returned.
    
    .. seealso::
    Save, 
    ExtractHistogram, 
    Merge, 
    NbEventSelect, 
    TimeScaling, 
    TimeSelect.
    """ 
    if filename !=None:
        return _Time_events(filename)
    else:
        raise TypeError("invalid filename")





