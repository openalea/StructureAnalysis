#!/usr/bin/env python
#-*- coding: utf-8 -*-
"""

.. topic:: enums.py summary

    A module to setup enums related to sequence analysis

    :Code status: mature
    :Documentation status: to be completed
    :Author: Thomas Cokelaer <Thomas.Cokelaer@sophia.inria.fr>

    :Revision: $Id: enums.py 15827 2014-02-19 16:17:45Z jbdurand $
"""
__version__ = "$Id: enums.py 15827 2014-02-19 16:17:45Z jbdurand $"

#from openalea.stat_tool._stat_tool import GraphicalType
from openalea.sequence_analysis._sequence_analysis import *

from openalea.stat_tool import  _stat_tool
#import openalea.sequence_analysis._sequence_analysis as _sequence_analysis

import openalea.stat_tool.enums as enumerate_st



#todo: get rid if type_dict, rename the enumerate properly

type_dict = enumerate_st.pearson_type



norm_type = {"Approximated": NormType.APPROXIMATED,
             "Exact": NormType.EXACT }

# not complete !! see MarkovianSequnece
markovian_sequence_type = {
    "SelfTransition": MarkovianSequenceType.SELF_TRANSITION,
    "Observation": MarkovianSequenceType.OBSERVATION,
    "Counting": MarkovianSequenceType.COUNTING,
    "Intensity": MarkovianSequenceType.INTENSITY,
    "Recurrence": MarkovianSequenceType.RECURRENCE_TIME,
    "Sojourn": MarkovianSequenceType.SOJOURN_TIME,
    "InitialRun": MarkovianSequenceType.INITIAL_RUN,
    "FinalRun": MarkovianSequenceType.FINAL_RUN,
    "Length": MarkovianSequenceType.LENGTH,
    "NbRun": MarkovianSequenceType.NB_RUN,
    "NbOccurrence": MarkovianSequenceType.NB_OCCURRENCE,
    "FirstOccurrence": \
        MarkovianSequenceType.FIRST_OCCURRENCE,
    "Mean": MarkovianSequenceType.SEQUENCE_MEAN,
    "Cumul": MarkovianSequenceType.SEQUENCE_CUMUL
        }


#todo: get rid if mode_type, rename the enumerate properly
mode_type = enumerate_st.round_type


#constants: 
NB_STATE_SEQUENCE = NB_STATE_SEQUENCE
NB_SEGMENTATION = NB_SEGMENTATION


INTER_EVENT = RenewalType.INTER_EVENT
WITHIN_OBSERVATION_PERIOD =  \
    RenewalType.WITHIN_OBSERVATION_PERIOD
LENGTH_BIAS =  RenewalType.LENGTH_BIAS
BACKWARD_RECURRENCE_TIME = \
    RenewalType.BACKWARD_RECURRENCE_TIME
FORWARD_RECURRENCE_TIME = RenewalType.FORWARD_RECURRENCE_TIME
NB_EVENT = RenewalType.NB_EVENT
MIXTURE = RenewalType.MIXTURE

renewal_nb_event_map = {
    "InterEvent":   RenewalType.INTER_EVENT,
    "LengthBias":   RenewalType.LENGTH_BIAS,
    "Backward":     RenewalType.BACKWARD_RECURRENCE_TIME,
    "Forward":      RenewalType.FORWARD_RECURRENCE_TIME,
    "Mixture":      RenewalType.MIXTURE,
    "Within":       RenewalType.WITHIN_OBSERVATION_PERIOD,
    "NbEvent":      RenewalType.NB_EVENT,
                        }

sub_func_map = {
    "Sequence": OutputType.SEQUENCE,
    "SubtractionResidual": OutputType.SUBTRACTION_RESIDUAL ,
    "DivisionResidual":OutputType.DIVISION_RESIDUAL,
    "Residual": OutputType.SUBTRACTION_RESIDUAL,

    }
func_map = {
    "Sequence": OutputType.SEQUENCE,
    "Trend": OutputType.TREND,
    "SubtractionResidual": OutputType.SUBTRACTION_RESIDUAL ,
    "Residual": OutputType.SUBTRACTION_RESIDUAL,
    "DivisionResidual":OutputType.DIVISION_RESIDUAL,
    "StandardizedResidual": OutputType.STANDARDIZED_RESIDUAL
    }
#todo: simplify
output_sequence_map = func_map
output_map = func_map

estimator_map = {
    "MaximumLikelihood": Estimator.MAXIMUM_LIKELIHOOD,
    "Laplace": Estimator.LAPLACE,
    "AdaptativeLaplace": Estimator.ADAPTATIVE_LAPLACE,
    "UniformSubset": Estimator.UNIFORM_SUBSET,
    "UniformCardinality": Estimator.UNIFORM_CARDINALITY
    }
estimator = estimator_map


SEQUENCE = 0
SUBTRACTION_RESIDUAL = 1
STANDARDIZED_RESIDUAL = 2

#todo: replace by appropriate enumerate from_se uence_analysis


seq_map = {
            "Observation":      _stat_tool.GraphicalType.OBSERVATION,
            "FirstOccurrence":  _stat_tool.GraphicalType.FIRST_OCCURRENCE,
            "Recurrence":       _stat_tool.GraphicalType.RECURRENCE_TIME,
            "Sojourn":          _stat_tool.GraphicalType.SOJOURN_TIME,
            "NbRun":            _stat_tool.GraphicalType.NB_RUN,
            "NbOccurrence":     _stat_tool.GraphicalType.NB_OCCURRENCE,
            "Forward": -1}

_INITIAL_RUN = 4
_FINAL_RUN = 5
_LENGTH = 8
_SEQUENCE_CUMUL = 9
_SEQUENCE_MEAN = 10

histogram_type = {
                  "FinalRun":5,
                  "InitialRun":4,
                  }

model_type = {
              "Multinomial": ChangeType.CATEGORICAL_CHANGE,
              "Poisson": ChangeType.POISSON_CHANGE,
              "Ordinal": ChangeType.ORDINAL_GAUSSIAN_CHANGE,
              "Gaussian": ChangeType.GAUSSIAN_CHANGE,
              "Mean": ChangeType.MEAN_CHANGE,
              "Variance": ChangeType.VARIANCE_CHANGE,
              "MeanVariance": ChangeType.MEAN_VARIANCE_CHANGE
              }

# todo:: ckean the 2 next enums ?
type_map = {
            "INT": INT_VALUE,
            "REAL" : REAL_VALUE,
            "STATE": STATE,
            # "NB_INTERNODE":  NB_INTERNODE,
            "AUXILIARY":  AUXILIARY,
            }

index_parameter_type_map = {
        "IMPLICIT_TYPE": IMPLICIT_TYPE,
        "TIME": TIME,
        "TIME_INTERVAL": TIME_INTERVAL,
        "POSITION": POSITION,
        "POSITION_INTERVAL": POSITION_INTERVAL
        }


stochastic_process_type = {
    'Ordinary': 'o',
    'Equilibrium' : 'e'
    }

markovian_algorithms = {
    'Forward':_stat_tool.RestorationAlgorithm.FORWARD,
    'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,
    'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,
    'ForwardBackwardSampling': \
        _stat_tool.RestorationAlgorithm.FORWARD_DYNAMIC_PROGRAMMING,
    'GeneralizedViterbi': \
        _stat_tool.RestorationAlgorithm.GENERALIZED_VITERBI,
    'Gibbs':_stat_tool.RestorationAlgorithm.GIBBS_SAMPLING,
    'NoComputation':_stat_tool.RestorationAlgorithm.NO_COMPUTATION,
    'Viterbi':_stat_tool.RestorationAlgorithm.VITERBI,
    }

sub_markovian_algorithms = {
    'EM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD,
    'MCEM':_stat_tool.RestorationAlgorithm.FORWARD_BACKWARD_SAMPLING,
    }


sub_markovian_algorithms_2 = {
    'Forward':_stat_tool.RestorationAlgorithm.FORWARD,
    'Viterbi':_stat_tool.RestorationAlgorithm.VITERBI,
    }


algorithm = {
    'CTM_BIC': Algorithm.CTM_BIC ,
    'BIC': Algorithm.CTM_BIC ,
    'CTM_KT': Algorithm.CTM_KT,
    'KT': Algorithm.CTM_KT,
    'LocalBIC': Algorithm.LOCAL_BIC,
    'Context': Algorithm.CONTEXT
}



likelihood_penalty_type = {
    'AIC': _stat_tool.LikelihoodPenaltyType.AIC,
    'AICc': _stat_tool.LikelihoodPenaltyType.AICc,
    'BIC': _stat_tool.LikelihoodPenaltyType.BIC,
    'BICc': _stat_tool.LikelihoodPenaltyType.BICc,
    'ICL' : _stat_tool.LikelihoodPenaltyType.ICL,
    'ICLc': _stat_tool.LikelihoodPenaltyType.ICLc,
    }


mean_computation_map = {
    "Computed" : _stat_tool.ClusterType.COMPUTED,
    "Estimated" : _stat_tool.ClusterType.ESTIMATED,
    "OneStepLate" : _stat_tool.ClusterType.ONE_STEP_LATE
    }

estimator_semi_markov_type = {
    "CompleteLikelihood" : _stat_tool.EstimatorHSMType.COMPLETE_LIKELIHOOD,
    "PartialLikelihood" :  _stat_tool.EstimatorHSMType.PARTIAL_LIKELIHOOD,
    "KaplanMeier" : _stat_tool.EstimatorHSMType.KAPLAN_MEIER
    }


ident_map = {
             "VOID" : _stat_tool.I_DEFAULT,
             "LINEAR": _stat_tool.RegressionType.STAT_LINEAR,
             "MONOMOLECULAR": _stat_tool.RegressionType.STAT_MONOMOLECULAR,
             "LOGISTIC": _stat_tool.RegressionType.STAT_LOGISTIC,
             "NONPARAMETRIC": _stat_tool.RegressionType.STAT_NONPARAMETRIC,

             }

output_type = {
               "ChangePoint" : ChangePointType.CHANGE_POINT,
               "Segment" : ChangePointType.SEGMENT
               }

output_display = {
               "ChangePoint" : ChangePointType.CHANGE_POINT,
               "Segment" : ChangePointType.SEGMENT,
               "State" : SemiMarkovState.SSTATE,
               "InState" : SemiMarkovState.IN_STATE,
               "OutState" : SemiMarkovState.OUT_STATE
            }


nb_segment_map = {
                  "Fixed": False,
                  "Estimated":True
                  }


begin_aligned_map = {
                  "Aligned": False,
                  "Free":True
                  }
end_aligned_map = begin_aligned_map

from openalea.sequence_analysis._sequence_analysis import \
    _Sequences, _MarkovianSequences, _VariableOrderMarkovData, \
    _SemiMarkovData, _NonHomogeneousMarkovData, _VariableOrderMarkov,\
    _HiddenVariableOrderMarkov, _HiddenSemiMarkov, _SemiMarkov, \
    _Renewal,_RenewalData,_TimeEvents, \
    _Correlation, _NonHomogeneousMarkovData
    # Tops, _TopParameters, \
    
    

sequence_alignment_first_arg = [_Sequences,
                                _MarkovianSequences,
                                _VariableOrderMarkovData,
                                _SemiMarkovData,
                                _NonHomogeneousMarkovData]

markov_model_comparison_first_arg = \
    [_VariableOrderMarkov,
     _HiddenVariableOrderMarkov,
     _HiddenSemiMarkov,
     _SemiMarkov]

markov_model_for_sequences_first_arg = \
    [_MarkovianSequences,
     _VariableOrderMarkovData,
     _SemiMarkovData,
     _NonHomogeneousMarkovData]

markov_model_for_sequences_second_arg = \
    [_VariableOrderMarkov,
     _SemiMarkov,
     _HiddenVariableOrderMarkov,
     _HiddenSemiMarkov]

ms_vomd_smd_list = [_MarkovianSequences,
                    _VariableOrderMarkovData,
                    _SemiMarkovData]

# todo: use those simpler naming convention:
ms_vomd_smd_nhmd = markov_model_for_sequences_first_arg

output_sequence = {
                   "DistanceMatrix": 'm',
                   "Sequences": 's'
                   }

indel_cost_map = {"Adaptative": IndelCost.ADAPTATIVE,
                     "Fixed": IndelCost.FIXED}


all_sequences_types = [ _VariableOrderMarkov,
                        _VariableOrderMarkovData,
                        _HiddenVariableOrderMarkov,
                        _HiddenSemiMarkov,
                        _SemiMarkov,
                        _MarkovianSequences,
                        _SemiMarkovData,
                        _NonHomogeneousMarkovData,
                        # _NonHomogeneousMarkov,
                        _Renewal,
                        _RenewalData,
                        _TimeEvents,
                        _Sequences,
                        # _Tops,
                        # _TopParameters,
                        _Correlation
                        ]

