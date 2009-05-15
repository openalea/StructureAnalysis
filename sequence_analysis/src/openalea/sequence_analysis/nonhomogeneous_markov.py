"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Nonhomogeneous_markov

import _sequence_analysis

__all__ = ['NonhomogeneousMarkov',
           '_Nonhomogeneous_markov']


# Extend dynamically class
interface.extend_class( _Nonhomogeneous_markov, interface.StatInterface)

# Add methods to _Vectors


def NonhomogeneousMarkov(*args, **kargs):
    """NonhomogeneousMarkov
   
    """ 

     
    Length = kargs.get("Length", 40)
    
    if isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            
            output = _Nonhomogeneous_markov(filename, Length)

        else:
            raise IOError("bad file name")
        
        
    return output





