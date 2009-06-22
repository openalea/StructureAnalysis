"""Semi markov

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/
"""
__revision__ = "$Id:  $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Semi_markov


__all__ = ['SemiMarkov',
           '_Semi_markov']


# Extend dynamically class
interface.extend_class( _Semi_markov, interface.StatInterface)

# Add methods to _Vectors


def SemiMarkov(filename=None, length=40, counting=True, cumul_threshold=0):
    """SemiMarkov
    
    Construction of a semi-Markov chain from an ASCII file.
    
    :Usage:

    >>> SemiMarkov(file_name, Length=40, Counting=True)    
  
    :Arguments:
    
    * file_name (string).
    
    :Optional Arguments: 
    
    * Length (int): length of sequences for the computation of the intensity and 
      counting characteristic distributions (default value: 20),
    * Counting (bool): computation of counting characteristic distributions default value: True).
    
    :Returned Object:
    
    If the construction succeeds, an object of type semi-markov is returned, 
    otherwise no object is returned.    
  
    :Background:
    
    A semi-Markov chain is constructed from a first-order Markov chain 
    representing transition between distinct states and state occupancy 
    distributions associated to the non-absorbing states. The state occupancy
    distributions are defined as object of type distribution with the
    additional constraint that the minimum time spent in a given state
    is at least 1 (inf_bound >= 1).
    
    .. seealso::
    
        :class:`~openalea.stat_tool.output.Save`,
        :func:`~openalea.sequence_analysis.compare.Compare`,
        :func:`~openalea.sequence_analysis.simulate.Simulate`.
 
    """ 
    if filename==None:
        raise TypeError('No filename provided.')
    elif not os.path.isfile(filename):
        raise OSError("Invalid filename")
    else:
        return _Semi_markov(filename, length, counting, cumul_threshold)





