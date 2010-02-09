"""NonhomogeneousMarkov

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id$"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Nonhomogeneous_markov
from openalea.stat_tool import error
from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH


__all__ = ['NonhomogeneousMarkov',
           '_Nonhomogeneous_markov']


# Extend dynamically class
interface.extend_class( _Nonhomogeneous_markov, interface.StatInterface)

# Add methods to _Vectors


def NonhomogeneousMarkov(*args, **kargs):
    """NonhomogeneousMarkov
   
    """ 

    error.CheckArgumentsLength(args, 1, 1)
    filename = args[0]
    Length = kargs.get("Length", DEFAULT_LENGTH)
    
    error.CheckType([filename, Length], [str, int])
    
        
    if os.path.isfile(filename):
        output = _Nonhomogeneous_markov(filename, Length)
    else:
        raise IOError("bad file name")
        
    return output





