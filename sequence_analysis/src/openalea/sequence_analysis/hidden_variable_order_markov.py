"""HiddenVariableOrderMarkov

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id$"

import os

import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import \
    _Hidden_variable_order_markov

from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH
from openalea.sequence_analysis._sequence_analysis import OCCUPANCY_THRESHOLD

from openalea.stat_tool import error

__all__ = ['HiddenVariableOrderMarkov',
           '_Hidden_variable_order_markov']


# Extend dynamically class
interface.extend_class( _Hidden_variable_order_markov, interface.StatInterface)

# Add methods to _Vectors


def HiddenVariableOrderMarkov(filename=None, Length=DEFAULT_LENGTH,
                              CumulThreshold=OCCUPANCY_THRESHOLD):
    """HiddenVariableOrderMarkov
   
    .. todo:: documentation
    """ 
        
    error.CheckType([filename, Length, CumulThreshold], [str, int, float])

    if not os.path.isfile(filename):
        raise IOError("Invalid filename")
    else: 
        return _Hidden_variable_order_markov(filename, Length,  CumulThreshold)
        





