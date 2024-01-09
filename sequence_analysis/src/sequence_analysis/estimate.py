#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: estimate.py summary

    A module dedicated to Estimate functions

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: estimate.py 15827 2014-02-19 16:17:45Z jbdurand $
"""
__revision__ = "$Id: estimate.py 15827 2014-02-19 16:17:45Z jbdurand $"


from openalea.stat_tool import error

# constants
from openalea.stat_tool._stat_tool import  \
    ONE_STEP_LATE, \
    COMPUTED, \
    COMPLETE_LIKELIHOOD,\
    PARTIAL_LIKELIHOOD,\
    I_DEFAULT, \
    D_DEFAULT, \
    ORDER

#maps
from openalea.stat_tool.estimate import estimator_type
from openalea.stat_tool.estimate import outside_type
from openalea.stat_tool.estimate import smoothing_penalty_type

from openalea.sequence_analysis.enums import estimator_semi_markov_type
from openalea.sequence_analysis.enums import ident_map
from openalea.sequence_analysis.enums import mean_computation_map
#from openalea.sequence_analysis.enums import markovian_algorithms
from openalea.sequence_analysis.enums import mean_computation_map
from openalea.sequence_analysis.enums import sub_markovian_algorithms
from openalea.sequence_analysis.enums import algorithm
from openalea.sequence_analysis.enums import estimator
from openalea.sequence_analysis.enums import likelihood_penalty_type
from openalea.sequence_analysis.enums import stochastic_process_type

# structure class
from openalea.sequence_analysis._sequence_analysis import\
    _HiddenSemiMarkov,\
    _HiddenVariableOrderMarkov,\
    _VariableOrderMarkov,\
    _RenewalData,\
    _TimeEvents
# _Tops,\

from openalea.stat_tool._stat_tool import _DiscreteMixture
from openalea.stat_tool._stat_tool import _Compound
from openalea.stat_tool._stat_tool import _FrequencyDistribution
from openalea.stat_tool._stat_tool import _Convolution
from openalea.stat_tool._stat_tool import _DiscreteParametric
from openalea.stat_tool._stat_tool import _DiscreteParametricModel
from openalea.stat_tool._stat_tool import _Distribution

__all__ = ['Estimate']
# to be checked or improve. Maybe the c++ code could be more explicit, i.e.,
# e switched to elementary and so on.



def _estimate_non_homogeneous_markov(obj, *args, **kargs):
    """
    Estimate switch hidden_variable_order_markov
    """
    #to be finalised given that nb_state is reachable. need test and examples
    #nb_state = obj.get_marginal(0).nb_value
    # error.CheckType(args, 2+nb_state)
    CountingFlag = kargs.get("Counting", True)

    ident = [ident_map[x] for x in args]

    return obj.nonhomogeneous_markov_estimation(ident, CountingFlag)


def _estimate_hidden_variable_order_markov(obj, *args, **kargs):
    """
    Estimate switch hidden_variable_order_markov
    """
    from openalea.sequence_analysis._sequence_analysis import \
        MIN_NB_STATE_SEQUENCE, \
        MAX_NB_STATE_SEQUENCE, \
        NB_STATE_SEQUENCE_PARAMETER
    from openalea.stat_tool._stat_tool import \
        FORWARD_BACKWARD, \
        FORWARD_BACKWARD_SAMPLING

    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    CommonDispersion = kargs.get("CommonDispersion", False)
    NbIteration = kargs.get("NbIteration", 80)
    Counting = kargs.get("Counting", True)
    StateSequence = kargs.get("StateSequence", True)
    Parameter = kargs.get("Parameter", NB_STATE_SEQUENCE_PARAMETER)
    MinNbSequence = kargs.get("MinNbStateSequence", MIN_NB_STATE_SEQUENCE)
    MaxNbSequence = kargs.get("MaxNbStateSequence", MAX_NB_STATE_SEQUENCE)
    Algorithm = error.ParseKargs(kargs, "Algorithm", 'EM', \
                                 sub_markovian_algorithms)

    error.CheckType([CommonDispersion, Counting, GlobalInitialTransition, NbIteration,
                     MinNbSequence, MaxNbSequence, Parameter, StateSequence],
                     [bool, bool, bool, int, int, int, [int, float], bool])

    error.CheckType([args[0]], [_HiddenVariableOrderMarkov])

    # sanity check on arguments
    # this one can be check only when Chain will be public and exported in
    # export_variable_order_markov
    # if type == 'e' and kargs.get("GlobalInitialTransition")" raise Error
    if Algorithm != sub_markovian_algorithms["MCEM"]:
        options = ["Parameter", "MaxNbStateSequence", "MinNbStateSequence"]
        for option in options:
            if kargs.get(option):
                raise ValueError("If % is provided, Algorithm cannot be MCEM" %
                                 option)

    if Algorithm == FORWARD_BACKWARD:
        hmarkov = obj.hidden_variable_order_markov_estimation(
                args[0], GlobalInitialTransition, CommonDispersion,
                Counting, StateSequence, NbIteration)

    elif Algorithm == FORWARD_BACKWARD_SAMPLING:
        hmarkov = obj.hidden_variable_order_markov_stochastic_estimation(
                        args[0], GlobalInitialTransition, CommonDispersion,
                        MinNbSequence, MaxNbSequence, Parameter, Counting,
                        StateSequence, NbIteration)

    return hmarkov

def _estimate_renewal_count_data(obj, itype, **kargs):
    """
    Estimate switch renewal_count_data
    """
    Type = 'v'
    error.CheckType([obj, itype], [[_TimeEvents, _RenewalData], str])
    if isinstance(itype, str):
        if itype == "Ordinary":
            Type = 'o'
        elif itype == "Equilibrium":
            Type = 'e'
        else:
            raise AttributeError("type must be Ordinary or Equilibrium")
    else:
        raise AttributeError("type must be Ordinary or Equilibrium")


    Estimator = error.ParseKargs(kargs, "Estimator",
                                 'Likelihood', estimator_type)

    NbIteration = kargs.get("NbIteration", I_DEFAULT)
    error.CheckType([NbIteration], [int])

    InitialInterEvent = kargs.get("InitialInterEvent", None)
    error.CheckType([InitialInterEvent], [[type(None), _DiscreteParametricModel,
                                           _DiscreteMixture, _Convolution, _Compound]])

    EquilibriumEstimator = error.ParseKargs(kargs, "EquilibriumEstimator",
                            'CompleteLikelihood', estimator_semi_markov_type)

    InterEventMean = error.ParseKargs(kargs, "InterEventMean",
                            'Computed', mean_computation_map)

    Penalty = error.ParseKargs(kargs, "Penalty", "SecondDifference",
                               smoothing_penalty_type)

    Outside = error.ParseKargs(kargs, "Outside", "Zero", outside_type)
    Weight = kargs.get("Weight", -1.)
    error.CheckType([Weight], [[int, float]])

    if Type != 'e':
        if kargs.get("EquilibriumEstimator"):
            raise Exception("EquilibriumEstimator cannot be used with type='e'")
        if kargs.get("InterEventMean"):
            raise Exception("InterEventMean be used with type='e'")

    if Estimator == estimator_type['PenalizedLikelihood']:
        if kargs.get("InterEventMean") is None:
            InterEventMean = ONE_STEP_LATE
        elif InterEventMean == COMPUTED:
            raise ValueError("""
                Incompatible options Estimator and InterEventMean""")
    else:
        if kargs.get("Penalty"):
            raise ValueError("""Incompatible options Penalty with type o""")
        if kargs.get("Weight"):
            raise ValueError("""Incompatible options Weight with type o""")
        if kargs.get("Outside"):
            raise ValueError("""Incompatible options Outside with type o""")

    if InitialInterEvent:
        #cast from InitialInterEvent to Mixture, Compound should be done

        if isinstance(InitialInterEvent, _DiscreteParametricModel):
            InitialInterEvent = _DiscreteParametric(InitialInterEvent)
        else:
            InitialInterEvent = _Distribution(InitialInterEvent)
        renew = obj.estimation_inter_event_type(Type, InitialInterEvent,
                                           Estimator, NbIteration,
                                           EquilibriumEstimator,
                                           InterEventMean, Weight,
                                           Penalty, Outside)
    else:
        renew = obj.estimation_type(Type, Estimator, NbIteration,
                               EquilibriumEstimator, InterEventMean ,
                               Weight, Penalty, Outside)

    return renew


def _estimate_renewal_interval_data(obj, **kargs):
    """
    Estimate switch renewal_count_data

    .. todo:: to be completed and validated with tests

    see stat_func4 in aml
    """
    #only LIKELIHOOD and PENALIZED_LIKELIHOOD
    Estimator = error.ParseKargs(kargs, "Estimator",
                                 'Likelihood', estimator_type)


    NbIteration = kargs.get("NbIteration", I_DEFAULT)
    error.CheckType([NbIteration], [int])

    # distribution
    InitialInterEvent = kargs.get("InitialInterEvent", None)
    error.CheckType([InitialInterEvent], [[type(None), _DiscreteParametricModel,
                                           _DiscreteMixture, _Convolution, _Compound]])
    if isinstance(InitialInterEvent, _DiscreteParametricModel):
        InitialInterEvent = _DiscreteParametric(InitialInterEvent)
    else:
        InitialInterEvent = _Distribution(InitialInterEvent)
    #cast initialInterEvent to parametric ?
    Penalty = error.ParseKargs(kargs, "Penalty", "SecondDifference",
                               smoothing_penalty_type)
    Weight = kargs.get("Weight", D_DEFAULT)
    error.CheckType([Weight], [[int, float]])
    Outside = error.ParseKargs(kargs, "Outside", "Zero", outside_type)
    error.CheckType([Weight], [[int, float]])

    InterEventMean = error.ParseKargs(kargs, "InterEventMean",
                            'Computed', mean_computation_map)


    if Estimator == estimator_type['PenalizedLikelihood']:
        if kargs.get("InterEventMean") is None:
            InterEventMean = ONE_STEP_LATE
        elif InterEventMean == COMPUTED:
            raise ValueError("""
                Incompatible options Estimator and InterEventMean""")
    else:
        if kargs.get("Penalty"):
            raise ValueError("""Incompatible options Penalty with type o""")
        if kargs.get("Weight"):
            raise ValueError("""Incompatible options Weight with type o""")
        if kargs.get("Outside"):
            raise ValueError("""Incompatible options Outside with type o""")

    if isinstance(obj, _FrequencyDistribution):
        if InitialInterEvent:
            renew = obj.estimation_inter_event(InitialInterEvent,
                                           Estimator, NbIteration,
                                           InterEventMean, Weight,
                                           Penalty, Outside)
        else:
            renew = obj.estimation(Estimator, NbIteration,
                                   InterEventMean ,
                                   Weight, Penalty, Outside)
    else:
        if InitialInterEvent:
            renew = obj.estimation_inter_event(InitialInterEvent,
                                           Estimator, NbIteration,
                                           InterEventMean, Weight,
                                           Penalty, Outside)
        else:
            renew = obj.estimation(Estimator, NbIteration,
                                   InterEventMean ,
                                   Weight, Penalty, Outside)

    return renew

def _estimate_semi_markov(obj, *args, **kargs):

    Type = 'v'
    #error.CheckType([args[0]], [str])

    Type = error.CheckDictKeys(args[0], stochastic_process_type)

    NbIteration = kargs.get("NbIteration", I_DEFAULT)
    Counting = kargs.get("Counting", True)
    Estimator = error.ParseKargs(kargs, "Estimator", "CompleteLikelihood",
                                 estimator_semi_markov_type)

    OccupancyMean = error.ParseKargs(kargs, "OccupancyMean",
                                      'Computed', mean_computation_map)

    error.CheckType([Counting, NbIteration], [bool, int])

    if Type != 'e' or Estimator == PARTIAL_LIKELIHOOD:
        if kargs.get(NbIteration):
            raise ValueError("Forbidden options Estimate NbIteration")
        if kargs.get("OccupancyMean"):
            raise ValueError("Forbidden options Estimate OccupancyMean")

    return obj.semi_markov_estimation(Type , Estimator , Counting,
                                        NbIteration , OccupancyMean)



def _estimate_hidden_semi_markov(obj, *args, **kargs):
    """
    .. doctest::
        :options: +SKIP

        >>> hsmc21 = Estimate(seq21, "HIDDEN_SEMI-MARKOV", hsmc0)

    """
    from openalea.sequence_analysis._sequence_analysis import \
        MIN_NB_STATE_SEQUENCE, \
        MAX_NB_STATE_SEQUENCE, \
        NB_STATE_SEQUENCE_PARAMETER



    from openalea.stat_tool._stat_tool import \
        FORWARD_BACKWARD, \
        FORWARD_BACKWARD_SAMPLING, \
        KAPLAN_MEIER

#    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    CommonDispersion = kargs.get("CommonDispersion", False)
    NbIteration = kargs.get("NbIteration", I_DEFAULT)
    Counting = kargs.get("Counting", True)
    StateSequence = kargs.get("StateSequence", True)
    Parameter = kargs.get("Parameter", NB_STATE_SEQUENCE_PARAMETER)
    MinNbSequence = kargs.get("MinNbStateSequence", MIN_NB_STATE_SEQUENCE)
    MaxNbSequence = kargs.get("MaxNbStateSequence", MAX_NB_STATE_SEQUENCE)
    Algorithm = error.ParseKargs(kargs, "Algorithm", 'EM', \
                                 sub_markovian_algorithms)
    Estimator = error.ParseKargs(kargs, "Estimator", 'CompleteLikelihood',
                                estimator_semi_markov_type)
    InitialOccupancyMean = kargs.get("InitialOccupancyMean", D_DEFAULT)
    MeanComputation = error.ParseKargs(kargs, "OccupancyMean", 'Computed',
                                      mean_computation_map)
    RandomInitialization = kargs.get("RandomInitialization", False)

    error.CheckType([CommonDispersion, Counting, NbIteration,
                     MinNbSequence, MaxNbSequence, Parameter, StateSequence,
                     InitialOccupancyMean],
                     [bool, bool, int, int, int, [int, float], bool,
                     [float, int]])

    if Algorithm != sub_markovian_algorithms["MCEM"]:
        options = ["Parameter", "MaxNbStateSequence", "MinNbStateSequence"]
        for option in options:
            if kargs.get(option):
                raise ValueError(
                    "If % is provided, Algorithm cannot be MCEM" % option)
    if Algorithm != sub_markovian_algorithms["EM"]:
        if Estimator == KAPLAN_MEIER:
            raise ValueError(
                "Estimator= KaplanMeier and Algorithm = MCEM not possible")


    error.CheckType([args[0]], [[str, _HiddenSemiMarkov]])
    if isinstance(args[0], str):
        Type = 'v'

        error.CheckType([args[1]], [int])
        NbState = args[1]

        if args[0] == "Ordinary":
            error.CheckArgumentsLength(args, 3, 3)
            error.CheckType([args[2]], [str])
            Type = 'o'
            if args[2] not in ["LeftRight", "Irreducible"]:
                raise ValueError(
                        "third argument must be LeftRight or Irreducible.")
            if args[2] == "LeftRight":
                LeftRight = True
            else:
                LeftRight = False
        elif args[0] == "Equilibrium":
            error.CheckArgumentsLength(args, 2, 2)
            Type = 'e'
            LeftRight = False
        else:
            raise AttributeError("type must be Ordinary or Equilibrium")

        if ((Type != 'e') or (Estimator == PARTIAL_LIKELIHOOD) or \
            (Algorithm != FORWARD_BACKWARD)) and \
            kargs.get(InitialOccupancyMean):
            raise ValueError("Incompatible user arguments")

        if Algorithm == FORWARD_BACKWARD:
            hsmarkov = obj.hidden_semi_markov_estimation_model( Type, NbState,
                         LeftRight, InitialOccupancyMean, CommonDispersion, Estimator,
                         Counting, StateSequence, NbIteration, MeanComputation, RandomInitialization)
            return hsmarkov

        elif Algorithm == FORWARD_BACKWARD_SAMPLING:
            hsmarkov = obj.hidden_semi_markov_stochastic_estimation_model(
                Type, NbState, LeftRight, InitialOccupancyMean, CommonDispersion,
                MinNbSequence, MaxNbSequence, Parameter, Estimator, Counting,
                StateSequence, NbIteration)
            return hsmarkov

    elif isinstance(args[0], _HiddenSemiMarkov):

        #todo: add these lines once Chain is public
        #if ((( (args[0].type == 'o')) or
        #     (Estimator == PARTIAL_LIKELIHOOD) or
        # (Algorithm != FORWARD_BACKWARD)) and \
        # kargs.get("InitialOccupancyMean")):
        #    raise ValueError("Incompatible arguments")

        hsmarkov = args[0]
        if Algorithm == FORWARD_BACKWARD:
            output = obj.hidden_semi_markov_estimation(hsmarkov,
                                CommonDispersion, Estimator, Counting,
                                StateSequence, NbIteration, MeanComputation)
            return output
        elif Algorithm == FORWARD_BACKWARD_SAMPLING:
            return obj.hidden_semi_markov_stochastic_estimation(hsmarkov,
                            CommonDispersion, MinNbSequence, MaxNbSequence,
                            Parameter, Estimator, Counting,
                            StateSequence, NbIteration)




def _estimate_variable_order_markov(obj, *args, **kargs):
    """
    EStimate on variable order markov

    """
    from openalea.sequence_analysis._sequence_analysis import \
        LOCAL_BIC_THRESHOLD,\
        CTM_KT_THRESHOLD,\
        CTM_BIC_THRESHOLD,\
        CONTEXT_THRESHOLD,\
        CTM_BIC,\
        CTM_KT,\
        CONTEXT,\
        LOCAL_BIC

    Order = kargs.get("Order", None)
    MaxOrder = kargs.get("MaxOrder", ORDER)
    MinOrder = kargs.get("MinOrder", 0)
    Threshold = kargs.get("Threshold", LOCAL_BIC_THRESHOLD)

    error.CheckType([Threshold, MaxOrder, MinOrder], [[int, float], int, int])


    Algorithm = error.ParseKargs(kargs, "Algorithm", "LocalBIC", algorithm)
    Estimator = error.ParseKargs(kargs, "Estimator", "Laplace", estimator)
    Penalty = error.ParseKargs(kargs, "Penalty", "BIC", likelihood_penalty_type)

    GlobalInitialTransition = kargs.get("GlobalInitialTransition", True)
    GlobalSample = kargs.get("GlobalSample", True)
    Counting = kargs.get("Counting", True)

    error.CheckType([Counting, GlobalSample, GlobalInitialTransition],
                    [bool, bool, bool])

    #args0 is a string
    if len(args)>0 and isinstance(args[0], str):
        Type = 'v'
        Type = error.CheckDictKeys(args[0], stochastic_process_type)

        # check validity of the input arguments following AML's code
        if Algorithm != LOCAL_BIC and not kargs.get("Threshold"):
            if Algorithm == CTM_BIC:
                Threshold = CTM_BIC_THRESHOLD
            elif Algorithm == CTM_KT:
                Threshold = CTM_KT_THRESHOLD
            elif Algorithm == CONTEXT:
                Threshold = CONTEXT_THRESHOLD
        if Algorithm == CTM_KT and kargs.get("Estimator"):
            raise ValueError("Forbidden combinaison of Algorithm and Estimator")

        order_estimation = True

        if Order is not None:
            order_estimation = False
            MaxOrder = Order

        if not order_estimation:
            options = ["Algorithm", "Estimator", "GlobalSample", "MinOrder",
                       "Threshold"]
            for option in options:
                if kargs.get(option):
                    raise ValueError("Order and %s cannot be used together" %
                                     option)
        if Type == 'e' and kargs.get("GlobalInitialTransition"):
            raise ValueError("""
            Type e and GlobalInitialTransition cannot be used together""")

        if order_estimation is True:
            markov = obj.variable_order_markov_estimation1(
                Type, MinOrder, MaxOrder, Algorithm, Threshold, Estimator ,
                  GlobalInitialTransition , GlobalSample , Counting)
        else:
            markov = obj.variable_order_markov_estimation2(
                    Type, MaxOrder, GlobalInitialTransition, Counting)

    #Variable order markov case
    elif isinstance(args[0], _VariableOrderMarkov):
        vom = args[0]
        # can be implemted once Chain class is public and exported
        # in export_variable_order_markov
     #   if vom.type == 'e' and kargs.get("GlobalInitialTransition"):
     #       raise ValueError("""
     #       Type e and GlobalInitialTransition cannot be used together""")

        markov = obj.variable_order_markov_estimation3(vom,
                      GlobalInitialTransition, Counting)

    # array case
    elif isinstance(args[0], list):
        symbol = args[0]
        markov = obj.lumpability_estimation(symbol, Penalty,
                                         Order, Counting)

    else:
        raise KeyError("jfjf")

    return markov


#todo: should be estimate_top_parameters ?
#def _estimate_top(obj, **kargs):
    #"""
    #Top parameters Estimate switch
    #"""
    #error.CheckType([obj], [_Tops])

    #MinPosition = kargs.get("MinPosition", 1)
    #MaxPosition = kargs.get("MaxPosition", obj.max_position)
    #Neighbourhood = kargs.get("Neighbourhood", 1)
    ## user may use us of uk spelling.
    #Neighbourhood = kargs.get("Neighborhood", Neighbourhood)
    #EqualProbability = kargs.get("EqualProbability", False)


    #error.CheckType([MinPosition, MaxPosition], [int, int])
    #return obj.estimation(MinPosition, MaxPosition, Neighbourhood,
                          #EqualProbability)



def _estimate_dispatch(obj, *args, **kargs):
    """


        fct_map = {
        "VARIABLE_ORDER_MARKOV": "estimate_variable_order_markov",
        "HIDDEN_VARIABLE_ORDER_MARKOV": "estimate_hidden_variable_order_markov",
        "HIDDEN_SEMI-MARKOV": "estimate_hidden_semi_markov"
        }
    """
    #fct = getattr(obj, "estimate_%s" % fct_map[type] )
    fct_map_distribution = [
        "NONPARAMETRIC" ,
        "NP" ,
        "B" ,
        "BINOMIAL" ,
        "P" ,
        "POISSON" ,
        "NB" ,
        "NEGATIVE_BINOMIAL" ,
        "U" ,
        "UNIFORM" ,
        "MIXTURE" ,
        "CONVOLUTION" ,
        "COMPOUND",
        "Equilibrium",
        "Ordinary",
        ]
    fct_map = ["VARIABLE_ORDER_MARKOV",
               "HIDDEN_VARIABLE_ORDER_MARKOV",
               "HIDDEN_SEMI-MARKOV",
               "SEMI-MARKOV",
               "NONHOMOGENEOUS_MARKOV",
               "MARKOV"
               ]

    if len(args) == 0:
        itype = None
    else:
        itype = args[0]

    if (itype not in fct_map_distribution) and (itype not in fct_map) and \
        (type(obj) not in [_TimeEvents, _RenewalData, _FrequencyDistribution]):
        raise KeyError("Valid type are %s or %s"
                       % (str(fct_map),
                          str(fct_map_distribution)))
    if itype == "VARIABLE_ORDER_MARKOV":
        return _estimate_variable_order_markov(obj, *args[1:], **kargs)
    elif itype == "HIDDEN_VARIABLE_ORDER_MARKOV":
        return _estimate_hidden_variable_order_markov(obj, *args[1:], **kargs)
    elif itype == "HIDDEN_SEMI-MARKOV":
        return _estimate_hidden_semi_markov(obj, *args[1:], **kargs)
    elif itype == "SEMI-MARKOV":
        return _estimate_semi_markov(obj, *args[1:], **kargs)
    elif itype == "NONHOMOGENEOUS_MARKOV":
        return _estimate_non_homogeneous_markov(obj, *args[1:], **kargs)
    elif (type(obj)==_TimeEvents) or (type(obj)==_RenewalData) \
            or (len(args)>=2 and isinstance(obj, _FrequencyDistribution) \
            and isinstance(args[0], _FrequencyDistribution)   \
            and isinstance(args[1], _FrequencyDistribution)):
        if len(args)>=1 and isinstance(obj, _FrequencyDistribution)==False and\
            type(args[0])==str:

            return  _estimate_renewal_count_data(obj, args[0], *args[1:], **kargs)
        else:

            #always 'equilibrium' as second argument
            return  _estimate_renewal_interval_data(obj, **kargs)
    else:
        from openalea.stat_tool.estimate import Estimate as HistoEstimate
        return HistoEstimate(obj, itype, *args[1:], **kargs)



def Estimate(obj, *args, **kargs):
    """Estimate

    * Estimation of distributions.
    * Estimation of 'top' parameters.
    * Estimation of a renewal process from count data.
    * Estimation of (hidden) Markovian models.

    :Usage:
    .. doctest::
        :options: +SKIP

        >>> Estimate(histo, "NON-PARAMETRIC")
        >>> Estimate(histo, "NB", MinInfBound=1, InfBoundStatus="Fixed")
        >>> Estimate(histo, "MIXTURE", "B", dist,..., MinInfBound=1,
            InfBoundStatus="Fixed", DistInfBoundStatus="Fixed")
        >>> Estimate(histo, "MIXTURE", "B", "NB",..., MinInfBound=1,
                InfBoundStatus="Fixed", DistInfBoundStatus="Fixed",
                NbComponent="Estimated", Penalty="AIC")
        >>> Estimate(histo, "CONVOLUTION", dist,MinInfBound=1, Parametric=False)
        >>> Estimate(histo, "CONVOLUTION", dist,InitialDistribution=initial_dist,
                Parametric=False)
        >>> Estimate(histo, "COMPOUND", dist, unknown, Parametric=False,
                MinInfBound=0)
        >>> Estimate(histo, "COMPOUND", dist, unknown,
                InitialDistribution=initial_dist, Parametric=False)

        >>> Estimate(top, MinPosition=1, MaxPosition=5, Neighbourhood=2,
                EqualProba=True)

        >>> Estimate(timev, type, NbIteration=10,Parametric=True)
        >>> Estimate(timev, type, InitialInterEvent=initial_dist,
                NbIteration=10, Parametric=True)

        >>> Estimate(seq, "NONHOMOGENEOUS_MARKOV", MONOMOLECULAR, VOID,
                Counting=False)
        >>> Estimate(seq, "SEMI-MARKOV", Counting=False)
        >>> Estimate(seq, "HIDDEN_MARKOV", nb_state, structure,
                SelfTransition=0.9, NbIteration=10,
                StateSequences="Viterbi", Counting=False)
        >>> Estimate(seq, "HIDDEN_MARKOV", hmc, Algorithm="Viterbi",
                NbIteration=10, Order=2, Counting=False)
        >>> Estimate(seq, "HIDDEN_MARKOV", "NbState", min_nb_state,
                max_nb_state, Penalty="AIC", Order=2, Counting=False)
        >>> Estimate(seq, "HIDDEN_MARKOV", "NbState", hmc, state,
                max_nb_state, Penalty="AIC", SelfTransition=0.9, Counting=False)
        >>> Estimate(seq, "HIDDEN_SEMI-MARKOV", nb_state, structure,
                OccupancyMean=20, NbIteration=10, Estimator="PartialLikelihood",
                StateSequences="Viterbi", Counting=False)
        >>> Estimate(seq, "HIDDEN_SEMI-MARKOV", hsmc, Algorithm="Viterbi",
            NbIteration=10, Counting=False)

    :Arguments:

    * histo (histogram, mixture_data, convolution_data, compound_data),
    * dist (distribution, mixture, convolution, compound),
    * unknown (string): type of unknown distribution: "Sum" or "Elementary".
    * top (tops),
    * timev (time_events, renewal_data),
    * type (string): type or renewal process: "Ordinary" or "Equilibrium".
    * seq (discrete_sequences, markov_data, semi-markov_data),
    * states, ... (array(int)):  new states corresponding to a partition of the
      original state space,
    * hmc (hidden_markov),
    * structure (string): structural properties of  the underlying Markov chain:
      "Irreductible" or "LeftRight" (i.e. a succession of transient states and a
      final absorbing state),
    * nb_state (int): number of states with 2 <= nb_state <=  15,
    * min_nb_state (int): minimum number of states,
    * max_nb_state (int): maximum number of states with 2 <= min_nb_state < max_nb_state <= 15
      or (number of states of the initial hidden Markov chain hmc) < max_nb_state<= 15.
    * state (int): state to be duplicated,
    * hsmc (hidden_semi-markov).

    :Optional Arguments:

    **distribution case**

    * MinInfBound (int): lower bound to the range of possible values (0 - default
      value - or 1). This optional argument cannot be used in conjunction with the
      optional argument InitialDistribution.
    * InfBoundStatus (string): shifting or not of the distribution: "Free" (default
      value) or "Fixed". This optional argument cannot be used if the second mandatory
      argument giving the model type is "NON-PARAMETRIC" ("NP").
    * DistInfBoundStatus (string): shifting or not of the subsequent components of
      the mixture: "Free" (default value) or "Fixed". This optional argument can
      only be used if the second mandatory argument giving the distribution type is "MIXTURE".
    * NbComponent (string): estimation of the number of components of the mixture:
      "Fixed" (default value) or "Estimated". This optional argument can only be
      used if the second mandatory argument giving the distribution type is "MIXTURE".
      the number of estimated components is comprised between 1 and a maximum number
      which is given by the number of specified parametric distributions in the
      mandatory arguments (all of these distributions are assumed to be unknown).
    * Penalty (string): type of penalty function for model selection: "AIC"
      (Akaike Information Criterion), "AICc" (corrected Akaike Information Criterion
      - default value) or "BIC" (Bayesian Information Criterion). This optional
      argument can only be used if the second mandatory argument giving the distribution
      type is "MIXTURE" and if the optional argument NbComponent is set at "Estimated".
    * Parametric (bool): reestimation of a discrete nonparametric or parametric
      distribution (default value: True). This optional argument can only be used if
      the second mandatory argument giving the distribution type is "CONVOLUTION"
      or "COMPOUND".
    * InitialDistribution (distribution, mixture, convolution, compound): initial
      distribution for the EM deconvolution-type algorithm. This optional argument
      can only be used if the second mandatory argument giving the distribution type
      is "CONVOLUTION" or "COMPOUND". This optional argument cannot be used in
      conjunction with the optional argument MinInfBound.

    .. note:: the optional arguments MinInfBound and InitialDistribution are mutually exclusive.


    **top case**

    * MinPosition (int): lower position taken into account for the estimation of 'top' parameters.
    * MaxPosition (int): upper position taken into account for the estimation of 'top' parameters.
    * Neighbourhood (int): neighbourhood taken into account for the estimation of 'top' parameters.
    * EqualProba (bool): growth probabilities of the parent shoot and of the offspring shoots equal or not (default value: False).

    **renewal case**

    * InitialInterEvent (distribution, mixture, convolution, compound): initial inter-event distribution for the EM algorithm.
    * NbIteration (int): number of iterations of the EM algorithm.
    * Parametric (bool): reestimation of a discrete nonparametric or parametric distribution (default value: False).

    **markovian case**

    * Counting (bool): computation of counting distributions (default value: True).
    * Order (int): Markov chain order (default value: 1). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV", "NONHOMOGENEOUS_MARKOV" or "HIDDEN_MARKOV".
    * MaxOrder (int): maximum order of the Markov chain (default value: 4). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV".
    * Penalty (string): type of penalty function for model selection: "AIC" (Akaike Information Criterion), "AICc" (corrected Akaike Information Criterion) or "BIC" (Bayesian Information Criterion). This optional argument can only be used if the second mandatory argument giving the model type is "MARKOV" (default value: "BIC") and if the optional argument MaxOrder is set or else, if a new set of states is given (defining a partition of the original state space) or else, if the second mandatory argument giving the model type is "HIDDEN_MARKOV" and the third "NbState" (default value: "AICc").
    * Algorithm (string): type of algorithm: "ForwardBackward" (the default) or "Viterbi". This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" or "HIDDEN_SEMI-MARKOV".
    * NbIteration (int): number of iterations of the estimation algorithm.
    * SelfTransition (real): self-transition probability. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" and if the initial model used in the iterative estimation procedure (EM algorithm) is only specified by its number of states, its structural properties and eventually its order.
    * OccupancyMean (int/real): average state occupancy. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_SEMI-MARKOV" and if the initial model used in the iterative estimation procedure (EM algorithm) is only specified by its number of states and its structural properties.
    * Estimator (string): type of estimator: "CompleteLikelihood" (the default) or "PartialLikelihood". In this latter case, the contribution of the time spent in the last visited state is not taken into account inthe estimation of the state occupancy distributions. This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_SEMI-MARKOV" and the optional argument Algorithm is set at "ForwardBackward".
    * StateSequences (string): Computation of the optimal state sequences: no computation (the default), "ForwardBackward" or "Viterbi". This optional argument can only be used if the second mandatory argument giving the model type is "HIDDEN_MARKOV" or "HIDDEN_SEMI-MARKOV" and if the optional argument Algorithm is not set at "Viterbi".


    :Returned Object:

    **distribution case**

    In case of success of the estimation procedure, the type of the returned object
    (chosen among distribution, mixture, convolution or compound) is given by the
    second mandatory argument. Otherwise no object is returned. The returned object
    of type distribution, mixture, convolution or compound contains both the estimated
    distribution and the data used in the estimation procedure. In the case of
    mixtures, convolutions, or compound (or stopped-sum) distributions, the returned
    object contains pseudo-data computed as a byproduct of the EM algorithm which
    can be extracted by the function ExtractData.

    **top case**

    In case of success of the estimation procedure, an object of type top_parameters
    is returned, otherwise no object is returned. The returned object of type top_parameters
    contains both the estimated model and the data used for the estimation.

    **renewal case**

    In case of success of the estimation procedure, an object of type renewal is
    returned, otherwise no object is returned. The returned object of type renewal
    contains both the estimated renewal process and the count data used in the
    estimation procedure.

    **markovian case**

    In case of success of the estimation procedure, the type of the returned object
    (chosen among markov, semi-markov, hidden_markov, hidden_semi-markov) is given
    by the second mandatory argument. Otherwise no object is returned. If the
    second mandatory argument is "NONHOMOGENEOUS_MARKOV", in case of success
    of the estimation procedure, the returned object is of type markov. If the
    second mandatory argument is "NONHOMOGENEOUS_MARKOV", the subsequent
    arguments chosen among "VOID" (homogeneous state), "MONOMOLECULAR" or
    "LOGISTIC", specify the evolution of the self-transition probabilities
    as a function of the index parameter. The returned object of type markov,
    semi-markov, hidden_markov or hidden_semi-markov contains both the estimated
    distribution and the data used in the estimation procedure. In the case of
    the estimation of a hidden Markov chain or a hidden semi-Markov chain,
    the returned object contains pseudo-data (optimal state sequences
    corresponding to the observed sequences used in the estimation procedure)
    computed as a byproduct of the EM algorithm which can be extracted by the
    function ExtractData.


    :Background:

    The aim of the model of 'tops' is to related the growth of offspring shoots
    to the growth of their parent shoot in the case of immediate branching. A
    model of 'tops' is defined by three parameters, namely the growth probability
    of the parent shoot, the growth probability of the offspring shoots (both in
    the sense of Bernoulli processes) and the growth rhythm ratio offspring
    shoots / parent shoot.

    :Description (markovian case):

    In the case of hidden Markovian models (second mandatory argument "HIDDEN_MARKOV"
    or "HIDDEN_SEMI-MARKOV"), either the forward-backward algorithm or the Viterbi
    algorithm can be used for estimation. The Viterbi algorithm should only be
    used for the estimation of hidden Markovian models based on an underlying
    "left-right" Markov chain (i.e. constituted of a succession of transient states
    and a final absorbing state). Hence, in this case, the model structure is
    implicitly "LeftRight" and should not be given as argument (only the number of
    states should be given as argument). Since the optimal state sequences are computed
    by the Viterbi algorithm, the optional argument StateSequences cannot be used if
    the optional argument Algorithm is set at "Viterbi".

    .. seealso::

        :func:`~openalea.stat_tool.data_transform.ExtractData`,
        :func:`~openalea.stat_tool.data_transform.ExtractDistribution`.
        :func:`~openalea.sequence_analysis.data_transform.AddAbsorbingRun`,
        :func:`ModelSelectionTest`.

    """

    # top case (no type specified and args may be empty)
    #if isinstance(obj, _Tops):
        #return _estimate_top(obj, *args, **kargs)
    #else:
    return _estimate_dispatch(obj, *args, **kargs)
