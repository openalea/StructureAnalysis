/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for StructureAnalysis developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#ifndef STAT_LABEL_H
#define STAT_LABEL_H

#if defined WIN32 || defined _WIN32 || defined __CYGWIN__
  #ifdef LIBSTAT_TOOL
    #ifdef __GNUC__
      #define EXTERN extern __attribute__ ((dllexport))
    #else
      #define EXTERN extern __declspec(dllexport)
    #endif
  #else
    #ifdef __GNUC__
      #define EXTERN extern __attribute__ ((dllimport))
    #else
      #define EXTERN extern __declspec(dllimport)
    #endif
  #endif
#else
  #if __GNUC__ >= 4
    #define EXTERN extern __attribute__ ((visibility ("default")))
  #else
    #define EXTERN extern 
  #endif
#endif

namespace stat_tool {


/****************************************************************
 *
 *  Key word identifiers (file format)
 */


  enum stat_tool_keyword {
    STATW_INF_BOUND ,
    STATW_SUP_BOUND ,
    STATW_PARAMETER ,
    STATW_PROBABILITY ,
    STATW_NO_SEGMENT ,
    STATW_SEQUENCE_LENGTH ,

    STATW_SHAPE ,
    STATW_SCALE ,
    STATW_ZERO_PROBABILITY ,
    STATW_MEAN ,
    STATW_STANDARD_DEVIATION ,
    STATW_MEAN_DIRECTION ,
    STATW_CONCENTRATION ,
    STATW_INTERCEPT ,
    STATW_SLOPE ,
    STATW_AUTOREGRESSIVE_COEFF ,

    STATW_MIXTURE ,
    STATW_CONVOLUTION ,
    STATW_DISTRIBUTIONS ,
    STATW_DISTRIBUTION ,
    STATW_WEIGHT ,
    STATW_WEIGHTS ,
    STATW_COMPONENTS ,
    STATW_COMPONENT ,

    STATW_COMPOUND ,
    STATW_SUM ,
    STATW_ELEMENTARY ,

    STATW_VARIABLE ,
    STATW_VARIABLES ,

    STATW_DISTANCE ,
    STATW_CATEGORIES ,
    STATW_PERIOD ,

    STATW_STATES ,
    STATW_INITIAL_PROBABILITIES ,
    STATW_TRANSITION_PROBABILITIES ,

    STATW_STATE ,
    STATW_FUNCTION ,

    STATW_OUTPUT_PROCESS ,
    STATW_OUTPUT_PROCESSES ,

    STATW_NONPARAMETRIC ,   // for ascending compatibility

    STATW_CATEGORICAL ,

    STATW_PARAMETRIC ,   // for ascending compatibility

    STATW_DISCRETE_PARAMETRIC ,
    STATW_CONTINUOUS_PARAMETRIC ,
    STATW_OBSERVATION_DISTRIBUTION ,
    STATW_OBSERVATION_MODEL ,

    STATW_OBSERVATION_PROBABILITIES ,  // for ascending compatibility

    STATW_OUTPUT
  };


  EXTERN const char *STAT_word[];
  EXTERN const char *STAT_discrete_distribution_word[];
  EXTERN const char *STAT_discrete_distribution_letter[];
  EXTERN const char *STAT_continuous_distribution_word[];
  EXTERN const char *STAT_continuous_distribution_letter[];
  EXTERN const char *STAT_variable_type_word[];
  EXTERN const char *STAT_variable_type_letter[];
  EXTERN const char *STAT_distance_word[];
  EXTERN const char *STAT_criterion_word[];
  EXTERN const char *STAT_function_word[];
  EXTERN const char *STAT_variable_word[];



/****************************************************************
 *
 *  Label identifiers
 */


  enum stat_tool_label {
    STATL_LINE ,
    STATL_WORD ,

    STATL_HIT_RETURN ,
    STATL_END ,

    STATL_MODEL ,

    STATL_MEAN ,
    STATL_MEDIAN ,
    STATL_MODE ,
    STATL_VARIANCE ,
    STATL_STANDARD_DEVIATION ,
    STATL_QUANTILE ,
    STATL_LOWER_QUARTILE ,
    STATL_UPPER_QUARTILE ,
    STATL_MEAN_ABSOLUTE_DEVIATION ,
    STATL_VARIATION_COEFF ,
    STATL_VARIANCE_MEAN_RATIO ,
    STATL_SKEWNESS_COEFF ,
    STATL_KURTOSIS_COEFF ,
    STATL_INFORMATION ,
    STATL_CONCENTRATION_COEFF ,
    STATL_SMOOTHNESS ,

    STATL_MEAN_DIRECTION ,
    STATL_MEAN_RESULTANT_LENGTH ,
    STATL_CIRCULAR_STANDARD_DEVIATION ,

    STATL_MEAN_CONFIDENCE_INTERVAL ,

    STATL_ONE_SIDED ,
    STATL_TWO_SIDED ,
    STATL_CHI2_TEST ,
    STATL_F_TEST ,
    STATL_T_TEST ,
    STATL_WILCOXON_MANN_WHITNEY_TEST ,
    STATL_KRUSKAL_WALLIS_TEST ,
    STATL_FREEDOM_DEGREE ,
    STATL_FREEDOM_DEGREES ,
    STATL_STANDARD_NORMAL_VALUE ,
    STATL_CHI2_VALUE ,
    STATL_F_VALUE ,
    STATL_T_VALUE ,
    STATL_CRITICAL_PROBABILITY ,
    STATL_REFERENCE ,
    STATL_TEST ,
    STATL_MANN_WHITNEY_INFERIOR_PROBABILITY ,
    STATL_MANN_WHITNEY_EQUAL_PROBABILITY ,
    STATL_MANN_WHITNEY_SUPERIOR_PROBABILITY ,

    STATL_FIT ,
    STATL_LIKELIHOOD ,
    STATL_NORMALIZED ,
    STATL_MAX_LIKELIHOOD ,
    STATL_DEVIANCE ,
    STATL_CLASSIFICATION_LIKELIHOOD ,
    STATL_MAX_CLASSIFICATION_LIKELIHOOD ,
    STATL_CLASSIFICATION_ENTROPY ,
    STATL_FREE_PARAMETER ,
    STATL_FREE_PARAMETERS ,
    STATL_PENALIZED_LIKELIHOOD ,
    STATL_ITERATION ,
    STATL_ITERATIONS ,

    STATL_DISTRIBUTION ,
    STATL_DISTRIBUTIONS ,
    STATL_UNPROPER ,
    STATL_COMPLEMENTARY_PROBABILITY ,
    STATL_CUMULATIVE ,
    STATL_SURVIVOR ,
    STATL_FUNCTION ,
    STATL_MATCHING ,
    STATL_CONCENTRATION ,
    STATL_CURVE ,

    STATL_MIXTURE ,
    STATL_WEIGHT ,
    STATL_COMPONENT ,
    STATL_CATEGORICAL ,
    STATL_DISCRETE_PARAMETRIC ,
    STATL_COMPONENT_DISTANCE ,
    STATL_POSTERIOR_ASSIGNMENT_PROBABILITY ,
    STATL_ASSIGNMENT_ENTROPY ,
    STATL_CONVOLUTION ,
    STATL_COMPOUND ,
    STATL_SUM ,
    STATL_ELEMENTARY ,
    STATL_DEATH_PROBABILITY ,
    STATL_SURVIVAL_PROBABILITY ,
    STATL_FREQUENCY ,

    STATL_BACKWARD ,
    STATL_FORWARD ,
    STATL_OBSERVATION_INTER_EVENT ,
    STATL_SOJOURN_TIME ,

    STATL_FREQUENCY_DISTRIBUTION ,
    STATL_FREQUENCY_DISTRIBUTIONS ,
    STATL_HISTOGRAM ,
    STATL_BIN_WIDTH ,
    STATL_SAMPLE ,
    STATL_SAMPLE_SIZE ,
    STATL_DISSIMILARITIES ,

    STATL_INFORMATION_RATIO ,
    STATL_CLUSTERING_STEP ,
    STATL_CATEGORY ,

    STATL_VECTOR ,
    STATL_VECTORS ,
    STATL_NB_VECTOR ,
    STATL_VARIABLE ,
    STATL_IDENTIFIER ,
    STATL_MIN_VALUE ,
    STATL_MAX_VALUE ,
    STATL_MARGINAL ,
    STATL_VARIANCE_COVARIANCE_MATRIX ,
    STATL_CORRELATION_MATRIX ,
    STATL_SPEARMAN_RANK_CORRELATION_MATRIX ,
    STATL_KENDALL_RANK_CORRELATION_MATRIX ,
    STATL_LIMIT_CORRELATION_COEFF ,
    STATL_SPEARMAN_LIMIT_RANK_CORRELATION_COEFF ,
    STATL_KENDALL_LIMIT_RANK_CORRELATION_COEFF ,
    STATL_SUP_NORM_DISTANCE ,
    STATL_OVERLAP ,

    STATL_CONTINGENCY_TABLE ,
    STATL_DEVIATION_TABLE ,
    STATL_CHI2_CONTRBUTION_TABLE ,

    STATL_VARIANCE_ANALYSIS ,
    STATL_VARIATION_SOURCE ,
    STATL_BETWEEN_SAMPLES ,
    STATL_WITHIN_SAMPLES ,
    STATL_TOTAL ,
    STATL_SQUARE_SUM ,
    STATL_MEAN_SQUARE ,

    STATL_INTERCEPT ,
    STATL_SLOPE ,
    STATL_EXPLANATORY_VARIABLE ,
    STATL_RESPONSE_VARIABLE ,
    STATL_CORRELATION_COEFF ,
    STATL_R_SQUARED ,
//    STATL_REGRESSION_VARIATION_TOTAL_VARIATION ,
    STATL_DETERMINATION_COEFF ,
    STATL_NULL_AUTOREGRESSIVE_COEFF_95_CONFIDENCE_LIMIT ,
    STATL_RESIDUAL ,
    STATL_STANDARDIZED_RESIDUAL ,
    STATL_LINEAR ,
    STATL_LINEAR_MODEL ,
    STATL_NONPARAMETRIC ,
    STATL_REGRESSION ,
    STATL_ESTIMATION ,

    STATL_DISTANCE ,
    STATL_CUMUL_DISTANCE ,
    STATL_LENGTH ,
    STATL_DELETION ,
    STATL_INSERTION ,
    STATL_INDEL ,
    STATL_MATCH ,
    STATL_SUBSTITUTION ,
    STATL_TRANSPOSITION ,
    STATL_RATE ,

    STATL_MATRIX ,
    STATL_NB_ROW ,
    STATL_NB_COLUMN ,

    STATL_SYMMETRY_RATE ,
    STATL_TRIANGLE_INEQUALITY_RATE ,

    STATL_CLUSTER ,
    STATL_CLUSTERS ,
    STATL_NEIGHBOR ,
    STATL_MOST_DISTANT ,
    STATL_DIAMETER ,
    STATL_SEPARATION ,
    STATL_PATTERN_LEVEL ,
    STATL_ISOLATED ,
    STATL_NON_ISOLATED ,
    STATL_WITHIN ,
    STATL_BETWEEN ,
    STATL_RATIO ,
    STATL_PROTOTYPE ,
    STATL_STEP ,
    STATL_CHILD ,
    STATL_COMPOSITION ,
    STATL_CHILD_CLUSTER_DISTANCE_COEFF ,
    STATL_DIAMETER_COEFF ,
    STATL_DENDROGRAM_SCALE ,

    STATL_CLASS ,
    STATL_STATE ,
    STATL_STATES ,
    STATL_TRANSIENT ,
    STATL_RECURRENT ,
    STATL_ABSORBING ,
    STATL_MEMORY ,
    STATL_STATIONARY_PROBABILITIES ,

    STATL_VALUE ,
    STATL_VALUES ,
    STATL_OUTPUT ,
    STATL_OUTPUT_PROCESS ,
    STATL_OBSERVATION ,
    STATL_OBSERVATION_PROBABILITIY_MATRIX ,
    STATL_THEORETICAL ,
    STATL_RESTORATION ,
    STATL_WEIGHTS ,
    STATL_OBSERVATION_DISTRIBUTION_DISTANCE ,
    STATL_CONSECUTIVE_STATE_OBSERVATION_DISTRIBUTION_DISTANCE ,
    STATL_Q_Q_PLOT
  };


  EXTERN const char *STAT_label[];



/****************************************************************
 *
 *  Identifiers of error messages for lexical analysis of files
 */


  enum stat_tool_parsing {
    STATP_KEYWORD ,
    STATP_FORMAT ,

    STATP_DISTRIBUTION_NAME ,
    STATP_PARAMETER_NAME ,
    STATP_PARAMETER_INDEX ,
    STATP_NB_PARAMETER ,
    STATP_SEPARATOR ,
    STATP_PARAMETER_VALUE ,
    STATP_PROBABILITY_SUM ,

    STATP_DATA_TYPE ,
    STATP_NB_TOKEN ,
    STATP_EMPTY_SAMPLE ,

    STATP_VALUE_ORDER ,
    STATP_MAX_VALUE ,

    STATP_NB_DISTRIBUTION ,
    STATP_DISTRIBUTION_INDEX ,
    STATP_WEIGHT_VALUE ,

    STATP_NB_COMPONENT ,
    STATP_COMPONENT_INDEX ,

    STATP_NB_VARIABLE ,
    STATP_VARIABLE_INDEX ,
    STATP_VARIABLE_TYPE ,

    STATP_NB_CATEGORY ,
    STATP_LOCAL_DISTANCE ,
    STATP_TRIANGLE_INEQUALITY ,
    STATP_PERIOD_VALUE ,

    STATP_NB_STATE ,
    STATP_ORDER ,
    STATP_INITIAL_PROBABILITY ,
    STATP_TRANSITION_PROBABILITY ,
    STATP_CHAIN_STRUCTURE ,
    STATP_IRREDUCIBLE ,

    STATP_STATE_INDEX ,
    STATP_OUTPUT_INDEX ,
    STATP_OBSERVATION_PROBABILITY ,
    STATP_NB_OBSERVATION_DISTRIBUTION ,
    STATP_NON_CONSECUTIVE_OUTPUTS ,
    STATP_OBSERVATION_DISTRIBUTION_OVERLAP ,
    STATP_OBSERVATION_DISTRIBUTION_NON_OVERLAP ,
    STATP_OBSERVATION_DISTRIBUTION_TYPE ,
    STATP_NB_OUTPUT_PROCESS ,
    STATP_OUTPUT_PROCESS_INDEX
  };


  EXTERN const char *STAT_parsing[];



/****************************************************************
 *
 *  Identifiers of error messages
 */


  enum stat_tool_error {
    STATR_FILE_NAME ,
    STATR_FILE_PREFIX ,
    STATR_EMPTY_SAMPLE ,
    STATR_SAMPLE_SIZE ,
    STATR_SAMPLE_SIZES ,
    STATR_NON_EXISTING_DISTRIBUTION ,
    STATR_NON_EXISTING_FREQUENCY_DISTRIBUTION ,
    STATR_NO_DATA ,

    STATR_PLOT_NB_DISTRIBUTION ,
    STATR_PLOT_NB_HISTOGRAM ,
    STATR_PLOT_NULL_VARIANCE ,

    STATR_VALUE_RANGE ,
    STATR_MIN_INF_BOUND ,
    STATR_ESTIMATION_FAILURE ,

    STATR_NB_COMPLETE_INTERVAL_TOO_SMALL ,
    STATR_COMPLETE_MIN_VALUE ,
    STATR_FORWARD_MIN_VALUE ,
    STATR_NO_EVENT_MIN_VALUE ,
    STATR_MEAN_ESTIMATION ,
    STATR_INTER_EVENT_SUPPORT ,

    STATR_DISTRIBUTION_INDEX ,
    STATR_FREQUENCY_DISTRIBUTION_INDEX ,
    STATR_NB_DISTRIBUTION ,
    STATR_WEIGHT_STEP ,
    STATR_DISTRIBUTION_FLAG ,

    STATR_NOT_PRESENT ,
    STATR_OUTPUT_PROCESS_INDEX ,
    STATR_OUTPUT_PROCESS_TYPE ,
    STATR_OPTIMAL_ASSIGNMENT ,

    STATR_VARIABLE_NB_VALUE ,
    STATR_NB_OUTPUT_PROCESS ,
    STATR_NB_OUTPUT ,
    STATR_POSITIVE_MIN_VALUE ,

    STATR_KNOWN_DISTRIBUTION_MIN_VALUE ,
    STATR_KNOWN_DISTRIBUTION_MEAN ,
    STATR_NB_ITERATION ,
    STATR_MIN_NB_ASSIGNMENT ,
    STATR_PENALTY_WEIGHT ,

    STATR_CLUSTERING_STEP ,
    STATR_NB_CLASS ,
    STATR_SHIFT_VALUE ,
    STATR_THRESHOLD_VALUE ,
    STATR_ROUNDED_VALUE ,
    STATR_SMALLER_THAN ,
    STATR_GREATER_THAN ,
    STATR_NOT_ALLOWED ,
    STATR_NB_CATEGORY ,
    STATR_NON_CONSECUTIVE_CATEGORIES ,
    STATR_MISSING_CATEGORY ,
    STATR_CLUSTER_LIMIT ,
    STATR_INFORMATION_RATIO ,
    STATR_NULL_INFORMATION ,
    STATR_SCALING_COEFF ,
    STATR_VALUE ,
    STATR_MIN_VALUE ,
    STATR_MAX_VALUE ,

    STATR_MARGINAL_FREQUENCY_DISTRIBUTION ,
    STATR_MARGINAL_HISTOGRAM ,
    STATR_HISTOGRAM_BIN_WIDTH ,
    STATR_HISTOGRAM_MIN_VALUE ,
    STATR_BAD ,
    STATR_NB_VECTOR ,
    STATR_NB_VARIABLE ,
    STATR_NB_SELECTED_VARIABLE ,
    STATR_NB_SUMMED_VARIABLE ,
    STATR_VARIABLE_TYPE ,
    STATR_VARIABLE_INDEX ,
    STATR_VARIABLE_INDICES ,
    STATR_ALREADY_USED ,
    STATR_ALREADY_SELECTED ,
    STATR_SAMPLE_INDEX ,
    STATR_RANK_CORRELATION_COMPUTATION ,
    STATR_SHIFTED_SCALED ,
    STATR_NB_VALUE_PERIOD ,
    STATR_NB_VALUE ,
    STATR_MISSING_VALUE ,
    STATR_LEAST_SQUARE_ALGORITHM ,
    STATR_REGRESSION_FAILURE ,

    STATR_MATRIX_DIMENSIONS ,
    STATR_INFINITE_DISTANCES ,
    STATR_SYMMETRICAL_MATRIX ,
    STATR_UNSYMMETRICAL_MATRIX ,
    STATR_UNNORMALIZED_DISSIMILARITY_MEASURES ,

    STATR_SYMMETRY ,
    STATR_TRIANGLE_INEQUALITY ,

    STATR_MATRIX_STRUCTURE ,
    STATR_SQUARE_MATRIX ,
    STATR_SINGLE_ELEMENT_CLUSTERS ,
    STATR_NB_CLUSTER ,
    STATR_PATTERN_TYPE ,
    STATR_PROTOTYPE_IDENTIFIER ,
    STATR_NUMBER ,

    STATR_NULL_INITIAL_PROBABILITY ,

    STATR_ODD ,
    STATR_NON_SYMMETRICAL_DISTRIBUTION ,
    STATR_UNPROPER_DISTRIBUTION ,
    STATR_NO_PERMUTATION
  };


  EXTERN const char *STAT_error[];


};  // namespace stat_tool



#endif
