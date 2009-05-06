"""Sequences"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.sequence_analysis._sequence_analysis import _Hidden_semi_markov

import _sequence_analysis

__all__ = ['HiddenSemiMarkov',
           '_Hidden_semi_markov']


# Extend dynamically class
interface.extend_class( _Hidden_semi_markov, interface.StatInterface)

# Add methods to _Vectors


def HiddenSemiMarkov(filename=None, length=40, counting=True, cumul_threshold=0):
    """HiddenSemiMarkov
   
    Construction of an object of type hidden_semi-markov from an ASCII file.
    
    :Usage:
    HiddenSemiMarkov(file_name, Length->40, Counting->False)    
   
    :Arguments:
    file_name (string)
    
    :Optional Arguments: 
    * Length (int): length of sequences for the computation of the intensity 
      and counting characteristic distributions (default value: 20),
    * Counting (bool): computation of counting characteristic distributions 
      (default value: True).
    
    :Returned Object:
    If the construction succeeds, an object of type hidden_semi-markov is 
    returned, otherwise no object is returned.    
  
    :Background:
    A hidden semi-Markov chain is constructed from an underlying semi-Markov
    chain (first-order Markov chain representing transition between distinct 
    states and state occupancy distributions associated to the nonabsorbing
    states) and nonparametric observation (or state-dependent) distributions.
    The state occupancy distributions are defined as object of type distribution 
    with the additional constraint that the minimum time spent in a given 
    state is 1 (inf_bound <= 1).
     
    
    .. seealso::
    Save, 
    Compare (Markovian models for sequences), 
    Compare (Markovian models), 
    Estimate (Markovian models), 
    ComputeStateSequences, 
    Simulate (Markovian models).
    """ 
#cumul_threshold=0
    if filename==None:
        raise TypeError('invalid filename')
    else:
        #todo: check old_format, cumulthreshold
        old_format = False
        cumul_threshold = 1 #???
        return _Hidden_semi_markov(filename, length, counting, cumul_threshold, 
                            old_format)





