"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Variable_order_markov
from openalea.sequence_analysis._sequence_analysis import _Sequences

import _sequence_analysis

__all__ = ['VariableOrderMarkov',
           '_Variable_order_markov']


# Extend dynamically class
interface.extend_class( _Variable_order_markov, interface.StatInterface)

# Add methods to _Vectors


def VariableOrderMarkov(*args, **kargs):
    """VariableOrderMarkov
    
    """ 
    
    DEFAULT_LENGTH = 20
    Length = kargs.get("Length", DEFAULT_LENGTH)
     
    if isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            vom =  _Variable_order_markov(filename, Length)
        else:
            raise IOError("bad file name")
        
        
    return vom




