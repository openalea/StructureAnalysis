"""Renewal

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id:  $"

import os
import openalea.stat_tool.interface as interface
from openalea.stat_tool._stat_tool import _Parametric
from openalea.sequence_analysis._sequence_analysis import _Renewal

__all__ = ['Renewal',
           '_Renewal']


# Extend dynamically class
interface.extend_class( _Renewal, interface.StatInterface)

# Add methods to _Vectors


def Renewal(*args, **kargs):
    """Renewal
    
    Construction of a (either ordinary or equilibrium) renewal process from an inter-event distribution or from an ASCII file.
    
    :Usage:

    >>> Renewal("BINOMIAL", inf_bound, sup_bound, proba,  Type="Equilibriun", ObservationTime=40)
    >>> Renewal("POISSON", inf_bound, param, Type="Equilibriun", ObservationTime=40)
    >>> Renewal("NEGATIVE_BINOMIAL", inf_bound, param, proba, Type="Equilibriun", ObservationTime=40) 
    >>> Renewal(inter_event, Type="Equilibriun", ObservationTime=40)
    >>> Renewal(file_name, Type="Equilibriun", ObservationTime=40)  
  
    :Arguments:

    * inf_bound (int): lower bound to the range of possible values (shift parameter),
    * sup_bound (int): upper bound to the range of possible values (only relevant for binomial or uniform distributions),
    * param (int, real): parameter of either the Poisson distribution or the negative binomial distribution.
    * proba (int, real): probability of 'success' (only relevant for binomial or negative binomial distributions).
    
    .. note:: the names of the parametric discrete distributions can be summarized by their first letters: "B" ("BINOMIAL"), "P" ("POISSON"), "NB" ("NEGATIVE_BINOMIAL").
    
    * inter_event (distribution, mixture, convolution, compound): inter-event distribution,
    * file_name (string).
    
    :Optional Arguments:
    
    * Type (string): type of renewal process: "Ordinary" or "Equilibriun" (the default).
    * ObservationTime (int): length of the observation period for the computation of the intensity and counting distributions (default value: 20),
    
    :Returned Object:

    If the construction succeeds, an object of type renewal is returned, otherwise no object is returned.    

    :Background:

    A renewal process is built from a discrete distribution termed the inter-event distribution which represents the time interval between consecutive events. Two types of renewal processes are available:
    * ordinary renewal process where the start of the observation period coincides with the occurrence time of an event (synchronism assumption),
    * equilibrium or stationary renewal process where the start of the observation period is independent of the process which generates the data (asynchronism assumption).
    In the case where the arguments are the name and the parameters of the inter-event distribution, the constraints on parameters described in the definition of the syntactic form of the type distribution apply (cf. File Syntax).
    
    .. seealso::
        :func:`~openalea.stat_tool.output.Save`, 
        :func:`~openalea.sequence_analysis.simulate.Simulate` (renewal process)
         
    .. todo :: ident should correspond to Binomail,B, NegativeBinomial and so on
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



