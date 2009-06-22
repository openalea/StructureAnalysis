"""HiddenVariableOrderMarkov

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id: r $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Hidden_variable_order_markov

import _sequence_analysis

__all__ = ['HiddenVariableOrderMarkov',
           '_Hidden_variable_order_markov']


# Extend dynamically class
interface.extend_class( _Hidden_variable_order_markov, interface.StatInterface)

# Add methods to _Vectors


def HiddenVariableOrderMarkov(Filename=None, Length=40,  CumulThreshold=0):
    """HiddenVariableOrderMarkov
   
   .. todo:: documentation
    """ 
#cumul_threshold=0
    if Filename == None:
        raise TypeError('invalid filename')
    else:
        #todo: check old_format, cumulthreshold
        CumulThreshold = 1 #???
        return _Hidden_variable_order_markov(Filename, Length,  CumulThreshold)





