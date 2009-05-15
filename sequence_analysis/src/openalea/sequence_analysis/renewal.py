"""Renewal"""
__revision__ = "$Id: vectors.py 6217 2009-04-08 12:40:15Z cokelaer $"

import os
import openalea.stat_tool.interface as interface
from openalea.stat_tool._stat_tool import _Parametric
from openalea.sequence_analysis._sequence_analysis import _Renewal

import _sequence_analysis

__all__ = ['Renewal',
           '_Renewal']


# Extend dynamically class
interface.extend_class( _Renewal, interface.StatInterface)

# Add methods to _Vectors


def Renewal(*args, **kargs):
    """
    .. todo:: to be tested. Any file example ? 
     
    ..todo :: ident should correspond to Binomail,B, NegativeBinomial and so on
    """ 
    
    type_map = {
        "Equilibirum":'e',
        "Ordinary": 'o'
    }
    
    Type = kargs.get("Type", "Equilibrium")
    ObservationTime = kargs.get("ObservationTime", 20)
    Scale = kargs.get("Scale", None)  #todo check default values !
        
    # a filename constructor. check that only one argument, which is a string
    # ------------------ todo ----------- not tested
    if len(args)==1 and isinstance(args[0], str):
        filename = args[0]
        if os.path.isfile(filename):
            renewal =  _Renewal(filename)
        else:
            raise IOError("bad file name")

    # otherwise, we switch to a cosntructor from a distribution
    elif isinstance(args[0], str):
    # ------------------ todo ----------- not tested
        inf_bound = 0
        sup_bound = 10
        parameter = 2
        probability = 0.5 
        RENEWAL_THRESHOLD = 1
        ident = 1 
        # ..todo :: ident should correspond to Binomail,B, NegativeBinomial and so on 
        inter_event = _Parametric(ident , inf_bound , sup_bound , parameter ,
                                   probability , RENEWAL_THRESHOLD);

        if Scale:
            scaled_inter_event = _Parametric(inter_event , Scale)
            renewal = Renewal(scaled_inter_event , Type , ObservationTime);
        else:
            renewal = _Renewal(inter_event , Type , ObservationTime)
      
 
    #    renewal = _Renewal(args[0], range(0,len(args[0])), index_parameter_type)
    # or may be provided by the user.
    else:
    # --------------------------  tested
        renewal = _Renewal(_Parametric(args[0]), Type, ObservationTime)
    
    
    
        
    return renewal 



