"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Time_events
from openalea.sequence_analysis._sequence_analysis import _Sequences

import _sequence_analysis

__all__ = ['TimeEvents',
           '_Time_events']


# Extend dynamically class
interface.extend_class( _Time_events, interface.StatInterface)

# Add methods to _Vectors


def TimeEvents(*args, **kargs):
    """TimeEvents
    
    Construction of data of type {time interval between two observation dates,
     number of events occurring between these two observation dates} from time 
     sequences, from an object of type HISTOGRAM or from an ASCII file.
     
    :Usage:
    
    >>> TimeEvents(seq1, begin_date, end_date, PreviousDate->3, NextDate->8)   
    >>> TimeEvents(seqn, variable, begin_date, end_date, PreviousDate->3, NextDate->8)   
    >>> TimeEvents(histo, time)
    >>> TimeEvents(file_name)
    
    >>> h = Histogram([1,1,1,2,2,2])
    >>> t = TimeEvents(h, 2)
        
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
    
    .. todo:: fix the build_time_events method to allows constructor with histogram
        issue: this method is in stat_tool and returns a time events so stat_tool requires to know sequence_analysis...
    """ 
    
    PreviousDate = kargs.get("PreviousDate", -1)
    NextDate = kargs.get("NextDate", -1)
     
    if len(args)==1 and isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            time_events =  _Time_events(filename)
        else:
            raise IOError("bad file name")
    elif isinstance(args[0], _Sequences):
        seq = args[0]
        nb_variable = seq.nb_variable()
        if nb_variable != 1:
            variable = args[0]
            begin_date = args[1]
            end_date = args[2]
        else:
            begin_date = args[0]
            end_date = args[1]
        
            
        time_events = seq.extract_time_events(variable, begin_date, end_date,
                                     PreviousDate, NextDate)

    else:
        #todo finish this code with examples ? 
        distribution = args[0]
        time = args[1]
        
        time_events = _Time_events(distribution, time)
        
             
        
        
    return time_events





def NbEventSelect(obj, min, max):

    return obj.nb_event_select(min, max)


