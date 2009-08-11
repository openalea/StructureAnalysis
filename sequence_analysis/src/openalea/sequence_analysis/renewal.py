"""Renewal

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
uthor:/

"""
__revision__ = "$Id:  $"

import os
import openalea.stat_tool.interface as interface
from openalea.stat_tool._stat_tool import _Parametric
from openalea.stat_tool._stat_tool import _ParametricModel
from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _Mixture

from openalea.sequence_analysis._sequence_analysis import _Renewal, _Renewal_data, _Sequences, _Markovian_sequences
from openalea.sequence_analysis._sequence_analysis import DEFAULT_TIME

from openalea.stat_tool import error
from openalea.stat_tool.enumerate import distribution_identifier_type

__all__ = ['Renewal',
           '_Renewal', '_Renewal_data', 'RenewalData']


# Extend dynamically class
interface.extend_class( _Renewal, interface.StatInterface)
interface.extend_class( _Renewal_data, interface.StatInterface)

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
    
    #todo: move this enym to enumerate.py
    type_map = {
        "Equilibrium":'e',
        "Ordinary": 'o'
    }
    
    Type = error.ParseKargs(kargs, "Type", "Equilibrium", type_map)
    ObservationTime = kargs.get("ObservationTime", DEFAULT_TIME)
    Scale = kargs.get("Scale", None)  #todo check default values !
    
    error.CheckType([args[0]], [[str, _Parametric]])    
    
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
    # ------------------ todo: ----------- not tested
        ident = distribution_identifier_type[args[0]]
        if args[0] == "BINOMIAL" or args[0] == "B":
            error.CheckArgumentsLength(args, 3, 3)
            inf_bound = args[1]
            sup_bound = args[2]
            probability = args[3]
            parameter = -1
        elif args[0] == "NEGATIVE_BINOMIAL" or args[0] == "NB":
            error.CheckArgumentsLength(args, 3, 3)
            inf_bound = args[1]
            sup_bound = -1
            parameter = args[2]
            probability = args[3] 
        elif args[0] == "Poisson" or args[0] == "P":
            error.CheckArgumentsLength(args, 2, 2)
            inf_bound = args[1]
            sup_bound = -1
            parameter = args[2]
            probability = args[3]
        else:
            raise NotImplemented("""case not implemented. First arg must be a 
                valid filename or a "BINOMIAL", "NEGATIVE_BINOMIAL, or "POISSON"
                """)  
             
        RENEWAL_THRESHOLD = 1
        ident = args[0] 
        #todo: check that this is the correct call. In aml call is parameter_input 
        inter_event = _Parametric(ident , inf_bound , sup_bound , parameter ,
                                   probability , RENEWAL_THRESHOLD);

        if Scale:
            scaled_inter_event = _Parametric(inter_event , Scale)
            renewal = Renewal(scaled_inter_event , Type , ObservationTime);
        else:
            renewal = _Renewal(inter_event , Type , ObservationTime)
      
 
    #    renewal = _Renewal(args[0], range(0,len(args[0])), index_parameter_type)
    # or may be provided by the user.
    elif type(args[0]) in [_Mixture, _Convolution, 
                           _Compound, _ParametricModel]:
    # --------------------------  tested
        renewal = _Renewal(_Parametric(args[0]), Type, ObservationTime)
    else:
         raise NotImplemented("check your first argument")
     
    return renewal 



def RenewalData(*args):

    error.CheckArgumentsLength(args, 2)
    
    error.CheckType([args[0]], [[_Sequences, _Markovian_sequences]])
    error.CheckType([args[1], args[2]], [int, int] )   
    
    
    seq = args[0]
    variable = seq.nb_variable
    begin_index = args[1]
    end_index = args[2]
    
    timev = seq.extract_renewal_data(variable, begin_index, end_index)
        
    return _Renewal_data(timev)
    