"""HiddenVariableOrderMarkov

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id: r $"


import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Hidden_variable_order_markov

from openalea.sequence_analysis._sequence_analysis import DEFAULT_LENGTH
from openalea.sequence_analysis._sequence_analysis import OCCUPANCY_THRESHOLD

from openalea.stat_tool import error

__all__ = ['HiddenVariableOrderMarkov',
           '_Hidden_variable_order_markov']


# Extend dynamically class
interface.extend_class( _Hidden_variable_order_markov, interface.StatInterface)

# Add methods to _Vectors


def HiddenVariableOrderMarkov(Filename=None, Length=DEFAULT_LENGTH,
                              CumulThreshold=OCCUPANCY_THRESHOLD):
    """HiddenVariableOrderMarkov
   
    .. todo:: documentation
    """ 
            
    if Filename:
        error.CheckType([Filename, Length, CumulThreshold], [str, int, float])
        return _Hidden_variable_order_markov(Filename, Length,  CumulThreshold)
    else:
        raise TypeError("Filename must be provided")
        





