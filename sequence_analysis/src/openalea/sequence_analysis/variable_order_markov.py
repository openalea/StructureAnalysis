"""Sequences"""
__revision__ = "$Id$"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _VariableOrderMarkov, _VariableOrderMarkovData
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




