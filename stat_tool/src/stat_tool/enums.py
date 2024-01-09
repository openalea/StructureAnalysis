"""common enumerates

:Author: Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
:Code status: mature
:Documentation status: to be completed
:Revision: $Id: enums.py 15183 2013-11-06 10:35:50Z jbdurand $
"""
__version__ = "$Id: enums.py 15183 2013-11-06 10:35:50Z jbdurand $"

from openalea.stat_tool._stat_tool import *

from openalea.stat_tool._stat_tool import _FrequencyDistribution, _DiscreteMixtureData, \
    _CompoundData,_ConvolutionData, _DiscreteDistributionData, _DiscreteMixture, \
    _Compound, _Convolution, _DiscreteParametricModel, _Distribution, \
    _Vectors, _Cluster, _DistanceMatrix, _VectorDistance, _Regression, \
    _MultivariateMixture, _MultivariateMixtureData


# map to enumerate in boost python
criterion_type = {
    "NearestNeighbor":  CriterionType.NEAREST_NEIGHBOR,
    "FarthestNeighbor": CriterionType.FARTHEST_NEIGHBOR,
    "Averaging":        CriterionType.AVERAGING,
    }

distance_type = {
    "ABSOLUTE_VALUE":   DistanceType.ABSOLUTE_VALUE,
    "QUADRATIC":        DistanceType.QUADRATIC
    }

fitting_type = {
    "STANDARD_NORMAL": FittingType.STANDARD_NORMAL,
    "CHI2": FittingType.CHI2,
    "FISHER": FittingType.FISHER,
    "STUDENT": FittingType.STUDENT
    }

round_type = {
    "FLOOR": RoundType.FLOOR,
    "ROUND": RoundType.ROUND,
    "CEIL": RoundType.CEIL
    }

distribution_identifier_type = {
    "B":                DistributionIdentifierType.BINOMIAL,
    "BINOMIAL":         DistributionIdentifierType.BINOMIAL,
    "P":                DistributionIdentifierType.POISSON,
    "POISSON":          DistributionIdentifierType.POISSON,
    "NB":               DistributionIdentifierType.NEGATIVE_BINOMIAL,
    "NEGATIVE_BINOMIAL": DistributionIdentifierType.NEGATIVE_BINOMIAL,
    "U": DistributionIdentifierType.UNIFORM,
    "UNIFORM": DistributionIdentifierType.UNIFORM,
    "M": DistributionIdentifierType.MULTINOMIAL,
    "MULTINOMIAL": DistributionIdentifierType.MULTINOMIAL,
    }

variable_type = {
    "O" : VariableType.ORDINAL,
    "ORDINAL" : VariableType.ORDINAL,
    "N" : VariableType.NUMERIC,
    "NUMERIC" : VariableType.NUMERIC,
    "S" : VariableType.SYMBOLIC,
    "SYMBOLIC" : VariableType.SYMBOLIC,
    "C" : VariableType.CIRCULAR,
    "CIRCULAR" : VariableType.CIRCULAR,
    }
sub_variable_type = {
    "N" : VariableType.NUMERIC,
    "NUMERIC" : VariableType.NUMERIC,
    "S" : VariableType.SYMBOLIC,
    "SYMBOLIC" : VariableType.SYMBOLIC,
    }

pearson_type = {
    "PEARSON": PearsonType.PEARSON,
    "SPEARMAN": PearsonType.SPEARMAN,
    "KENDALL": PearsonType.KENDALL,
    "SPEARMAN2": PearsonType.SPEARMAN2,
    "Pearson": PearsonType.PEARSON,
    "Spearman": PearsonType.SPEARMAN,
    "Kendall": PearsonType.KENDALL,
    "Spearman2": PearsonType.SPEARMAN2
    }

smoothing_penalty_type = {
    'FirstDifference': FIRST_DIFFERENCE,
    'SecondDifference': SECOND_DIFFERENCE,
    'Entropy': ENTROPY,
    }

estimator_type = {
    "Zero":                 EstimatorType.ZERO,
    "Likelihood":           EstimatorType.LIKELIHOOD,
    "PenalizedLikelihood":  EstimatorType.PENALIZED_LIKELIHOOD,
    "Parametric": EstimatorType.PARAMETRIC_REGULARIZATION,
    "ParametricRegularization": 
        EstimatorType.PARAMETRIC_REGULARIZATION
    }

outside_type = {
    "Zero": ZERO,
    "Continuation": CONTINUATION,
    }

#  enum_<stat_tool::wrap_util::UniqueInt<3, 11> >("TOBEDONEType")
#  .value("PARTIAL_LIKELIHOOD", PARTIAL_LIKELIHOOD)
#  .value("COMPLETE_LIKELIHOOD", COMPLETE_LIKELIHOOD)
#  .value("KAPLAN_MEIER", KAPLAN_MEIER)
#  .export_values()
#  ;


#  enum_<stat_tool::wrap_util::UniqueInt<3, 12> >("ClusterType")
#  .value("COMPUTED", COMPUTED)
#  .value("ESTIMATED", ESTIMATED)
#  .value("ONE_STEP_LATE", ONE_STEP_LATE)
#  .export_values()
#  ;

likelihood_penalty_type = {
    'AIC': AIC,
    'AICc': AICc,
    'BIC': BIC,
    'BICc': BICc,
    'ICL' : ICL,
    'ICLc': ICLc,
    }

algorithm_type = {
     "Agglomerative":   AlgoType.AGGLOMERATIVE,
     "Divisive":        AlgoType.DIVISIVE,
     "Ordering":        AlgoType.ORDERING,
     }




# used by Clustering
format_type = { \
    "ASCII" :'a',
    "SpreadSheet": 's',
    "" : 'n'
    }

cluster_type = \
    {
     "Step":        "cluster_step",
     "Limit":       "cluster_limit",
     "Information": "cluster_information"
     }

compound_type = {
    'Elementary': 'e',
    'Sum': 's',
    }

variance_type = {
    "N": NUMERIC,
    "NUMERIC" : NUMERIC,
    "O": ORDINAL,
    "ORDINAL":  ORDINAL,
    }

keep_type = {
             "keep": "Keep",
             "Keep": "Keep",
             "reject": "Reject",
             "Reject": "Reject",
             "Rejected": "Reject"
             }

bool_type = {
             "False": False,
             False:False,
             True:True,
             "True":True
             }


algo_map = {"Averaging": 'a',
            "LeastSquares": 's'
            }


histogram_types = [ _FrequencyDistribution, 
                    _DiscreteMixtureData, 
                    _CompoundData, 
                    _ConvolutionData,
                    _DiscreteDistributionData]


model_distribution_types = [ _DiscreteMixture,
                            _Compound, 
                            _Convolution, 
                            _DiscreteParametricModel, 
                            _Distribution]

all_stat_tool_types  = [_FrequencyDistribution,
                        _DiscreteMixtureData,
                        _CompoundData,
                        _ConvolutionData,
                        _DiscreteDistributionData,
                        _DiscreteMixture,
                        _MultivariateMixture,
                        _MultivariateMixtureData,
                        _Compound,
                        _Convolution,
                        _DiscreteParametricModel,
                        _Distribution,
                        _VectorDistance,
                        _Vectors,
                        _Cluster,
                        _DistanceMatrix,
                        _Regression
                        ]


output_display = {
               "ChangePoint" : None,
               "Segment" : None,
               "State" : None,
               "InState" : None,
               "OutState" : None
            }

