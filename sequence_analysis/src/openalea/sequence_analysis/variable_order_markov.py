"""Sequences"""
__revision__ = "$Id:  $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Variable_order_markov

from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH

from openalea.stat_tool import error

__all__ = ['VariableOrderMarkov',
           '_Variable_order_markov']


# Extend dynamically class
interface.extend_class( _Variable_order_markov, interface.StatInterface)

# Add methods to _Vectors


def VariableOrderMarkov(*args, **kargs):
    """VariableOrderMarkov
        
    """ 
    error.CheckArgumentsLength(args, 1, 1)
    error.CheckType([args[0]], [str])
    filename = args[0]
    Length = kargs.get("Length", DEFAULT_LENGTH) 
    
    if os.path.isfile(filename):
        vom =  _Variable_order_markov(filename, Length)
    else:
        raise IOError("bad file name %s" % filename)
        
        
    return vom




