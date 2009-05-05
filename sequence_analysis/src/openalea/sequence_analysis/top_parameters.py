"""Tops"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Top_parameters

__all__ = ['Top_parameters',
           '_Top_parameters']


# Extend dynamically class
interface.extend_class( _Top_parameters, interface.StatInterface)

# Add methods to _Vectors


def Top_parameters(*args, **kargs):
    """todo

    """
    
    if len(args)==1: 
        #filename case
        if isinstance(args[0], str):
            filename = args[0]
            if os.path.isfile(filename):
                return _Top_parameters(filename, True)
            else:
                raise IOError("bad file name")
        #sequences case
    else:
        raise TypeError("Expected a valid filename or a list of lists (e.g., [[1,0],[0,1]])")
    
    

 
    


