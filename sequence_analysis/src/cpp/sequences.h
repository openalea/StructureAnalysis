/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Exploring and Analyzing Plant Architecture
 *
 *       Copyright 1995-2018 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: sequences.h 18659 2015-11-03 07:28:02Z guedon $
 *
 *       Forum for V-Plants developers:
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



#ifndef SEQUENCES_H
#define SEQUENCES_H


#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
// #include "stat_tool/vectors.h"
#include "stat_tool/regression.h"


namespace sequence_analysis {



/****************************************************************
 *
 *  Constants
 */


  const int DEFAULT_LENGTH = 20;         // default sequence length

  const int SEQUENCE_NB_VARIABLE = 30;   // maximum number of variables

  const int PLOT_NB_SEQUENCE = 200;      // maximum number of plotted sequences (Gnuplot output)
  const int PLOT_LEGEND_NB_SEQUENCE = 15;  // maximum number of labeled sequences (Gnuplot output)
  const double GROWTH_FACTOR = 1.;       // growth factor for computing the first relative growth rate

  enum index_parameter_type {
    IMPLICIT_TYPE ,                      // implicit index parameter
    TIME ,                               // time
    TIME_INTERVAL ,                      // time interval
    POSITION ,                           // position
    POSITION_INTERVAL                    // between-position interval
  };

  enum sequence_pattern {
    LENGTH_PATTERN ,
    SEQUENCE_CUMUL ,
    SEQUENCE_MEAN ,
    FIRST_OCCURRENCE_PATTERN ,
    SOJOURN_TIME_PATTERN ,
    NB_RUN_PATTERN ,
    NB_OCCURRENCE_PATTERN
  };

  enum run_position {
    BEGIN_RUN ,
    END_RUN ,
  };

  enum correlation_variable_type {
    OBSERVED_VALUE ,
    OBSERVED_STATE ,
    THEORETICAL_STATE ,
    OBSERVED_OUTPUT ,
    THEORETICAL_OUTPUT
  };

  enum autocorrelation_function_type {
    AUTOREGRESSIVE ,
    WHITE_NOISE ,
    VOID
  };

  enum sequence_transformation {
    SEQUENCE_COPY ,
    REVERSE ,
    ADD_STATE_VARIABLE ,
    EXPLICIT_INDEX_PARAMETER ,
    REMOVE_INDEX_PARAMETER
  };

  enum categorical_sequence_process_transformation {
    CATEGORICAL_SEQUENCE_PROCESS_COPY ,
    INIT_OCCUPANCY
  };

  enum initial_run {
    UNCHANGED ,
    ADD_INITIAL_RUN ,
    REMOVE_INITIAL_RUN
  };

  enum memory_tree_selection {
    CTM_BIC ,                            // Context Tree Maximizing/BIC algorithm
    CTM_KT ,                             // Context Tree Maximizing/Krichevsky-Trofimov algorithm
    LOCAL_BIC ,                          // recursive pruning/BIC algorithm
    CONTEXT                              // Context algorithm
  };

  enum transition_estimator {
    MAXIMUM_LIKELIHOOD ,
    LAPLACE ,
    ADAPTATIVE_LAPLACE ,
    UNIFORM_SUBSET ,
    UNIFORM_CARDINALITY
  };

  enum segment_model {
    CATEGORICAL_CHANGE ,
    POISSON_CHANGE ,
    NEGATIVE_BINOMIAL_0_CHANGE ,
    NEGATIVE_BINOMIAL_1_CHANGE ,
    ORDINAL_GAUSSIAN_CHANGE ,
    GAUSSIAN_CHANGE ,
    MEAN_CHANGE ,
    VARIANCE_CHANGE ,
    LINEAR_MODEL_CHANGE ,
    INTERCEPT_SLOPE_CHANGE ,
    AUTOREGRESSIVE_MODEL_CHANGE ,
    STATIONARY_AUTOREGRESSIVE_MODEL_CHANGE ,
    BAYESIAN_POISSON_CHANGE ,
    BAYESIAN_GAUSSIAN_CHANGE
  };

  const double PRIOR_VARIANCE_FACTOR = 100.;  // factor for deducing the variance of
                                              // the gamma prior distribution
  const double PRIOR_SAMPLE_SIZE = 1.;  // hyperparameter of the Gaussian-gamma prior distribution
  const double PRIOR_DEGREES_OF_FREEDOM = 2.;  // hyperparameter of the Gaussian-gamma prior distribution
  const double PRIOR_DISPERSION_FACTOR = 10.;  // factor for deducing the sum of squared deviations
                                               // of the prior distribution
  const int PRIOR_PRECISION = 2;         // hyperparameter precision

  const double MAX_NB_WORD = 1.e7;       // maximum number of words

  const double STATIONARY_PROBABILITY_THRESHOLD = 1.e-8;  // threshold for the computation of the stationary distribution
                                                          // of an equilibrium model
  const int STATIONARY_PROBABILITY_LENGTH = 10000;  // maximum length for the computation of the stationary distribution
                                                    // of an equilibrium model
  const double LEAVE_INCREMENT = 1.e-6;  // threshold for stopping the computation of the probability
                                         // of leaving definitively a state/observed category

  const double CTM_BIC_THRESHOLD = 6.;   // threshold for memory pruning
  const double CTM_KT_THRESHOLD = 12.;   // threshold for memory pruning
  const double LOCAL_BIC_THRESHOLD = 10.;  // threshold for memory pruning
  const double CONTEXT_THRESHOLD = 5.;  // threshold for memory pruning

  const int MIN_NB_STATE_SEQUENCE = 1;   // number of state sequences, 1st iteration of the MCEM algorithm
  const int MAX_NB_STATE_SEQUENCE = 10;  // maximum number of state sequences (MCEM algorithm)
  const double NB_STATE_SEQUENCE_PARAMETER = 1.;  // parameter for the number of state sequences (MCEM algorithm)

  const double OCCUPANCY_THRESHOLD = 0.99999;  // threshold on the cumulative distribution function for determining
                                               // the upper bound of the support of a state occupancy distribution
  const double OCCUPANCY_MEAN = 10.;     // mean state occupancy

  const int POSTERIOR_PROBABILITY_NB_SEQUENCE = 300;  // maximum number of sequences for the output of
                                                      // the posterior probabilities of the most probable state sequences

  const int NB_STATE_SEQUENCE = 10;      // number of computed state sequences

  const int COUNTING_FREQUENCY_MAX_LENGTH = 10000;  // maximum sequence length for the extraction of the counting frequency distributions
  const int COUNTING_MAX_LENGTH = 500;   // maximum sequence length for the computation of the counting distributions

  const int NB_SEQUENCE = 100000;        // maximum number of generated sequences
  const int MAX_LENGTH = 1000000;        // maximum generated sequence length
  const int CUMUL_LENGTH = 1000000;      // maximum cumulative generated sequence length

  const double RESIDUAL_STANDARD_DEVIATION_COEFF = 1.e-6;  // threshold for the estimation of the residual standard deviation

  enum sequence_type {
    SEQUENCE ,
    SEQUENCE_SAMPLE ,
    TREND ,
    SUBTRACTION_RESIDUAL ,
    ABSOLUTE_RESIDUAL ,
    DIVISION_RESIDUAL ,
    STANDARDIZED_RESIDUAL ,
    SEGMENTATION_ENTROPY ,
    SEGMENTATION_DIVERGENCE ,
    LOG_LIKELIHOOD_SLOPE
  };

  enum output_sequence_format {
    LINE ,
    COLUMN ,
    VECTOR ,
    ARRAY ,
    POSTERIOR_PROBABILITY
  };

  enum insertion_deletion_cost {
    ADAPTATIVE ,
    FIXED
  };

  enum edit_operation { 
    DELETION ,
    BEGIN_END_DELETION ,
    INSERTION ,
    BEGIN_END_INSERTION ,
    MATCH ,
    SUBSTITUTION ,
    TRANSPOSITION
  };

  enum multiple_alignment {
    DATA ,
    GAP ,
    BEGIN_END_GAP
  };

  const int NB_ALIGNMENT = 1000000;      // maximum number of alignments
  const int DISPLAY_NB_ALIGNMENT = 30;   // maximum number of alignments for the detailed displayed output
  const int FILE_NB_ALIGNMENT = 300;     // maximum number of alignments for the detailed file output

  const double INDEL_FACTOR_1 = 0.51;    // factor for deducing the insertion/deletion cost - simple alignement
  const double INDEL_FACTOR_N = 0.51;    // factor for deducing le cout d'insertion/deletion cost - multiple alignement
  const double TRANSPOSITION_FACTOR = 0.;  // factor for deducing the transposition cost
  const double INDEL_DISTANCE = 1.0;     // insertion/deletion cost

  enum change_point_profile {
    CHANGE_POINT ,
    SEGMENT
  };

  const double CHANGE_POINT_UNCERTAINTY_PROBABILITY = 0.05;  // probability for defining change-point uncertainty intervals
  const int SEQUENCE_MAX_NB_COLUMN = 20; // threshold for the writing of sequences in text files
  const double ROUNDOFF_ERROR = 1.e-10;  // error on a sum of doubles
  const int PENALTY_SHAPE_SCALING_FACTOR = 100;  // scaling factor for the plot of log-likelihoods
  const int NB_SEGMENTATION = 10;        // number of computed segmentations
  const int SLOPE_NB_SEGMENT_RANGE = 5;  // minimum number of points for computing the log-likelihood slope
  const int DIMENSION_JUMP_NB_SEGMENT = SLOPE_NB_SEGMENT_RANGE + 2;  // minimum number of segments for the dimension jump method
  const double SLOPE_STEP = 0.01;        // slope step for the dimension jump method
  const double MAX_SLOPE = 100.;         // maximum slope for the dimension jump method
  const int MIN_DIMENSION_JUMP = 2;      // minimum dimension jump
  const double MIN_RANK_SQUARE_SUM = 1.e-2;  // default value in the case of a unique rank in a segment

  enum correlation_normalization {
    APPROXIMATED ,
    EXACT
  };

  const double CORRELATION_FREQUENCY_RATIO = 0.1; // frequency ratio for stopping the computation of
                                                  // the correlation function
  const double AUTOCORRELATION_FREQUENCY_RATIO = 0.2; // frequency ratio for stopping the computation of
                                                      // the state autocorrelation function
  const int AUTOCORRELATION_MIN_FREQUENCY = 20;  // minimum frequency for stopping the computation of
                                                 // the state autocorrelation function

  const int MAX_DIFFERENCING_ORDER = 3;  // maximum differentiation order
  const int POINTWISE_AVERAGE_NB_SEQUENCE = 250;  // maximum number of written sequences for the file output
  const int ABSORBING_RUN_LENGTH = 5;    // default length of the final absorbing run
  const int MAX_ABSORBING_RUN_LENGTH = 20;  // maximum length of the final absorbing run
  const double ABSORBING_RUN_STANDARD_DEVIATION_FACTOR = 10.;  // factor for defining the standard deviation of
                                                               // true final absorbing runs



/****************************************************************
 *
 *  Class definition
 */


  class VariableOrderMarkovChain;
  class VariableOrderMarkov;
  class VariableOrderMarkovData;
  class VariableOrderMarkovIterator;
  class HiddenVariableOrderMarkov;
  class SemiMarkov;
  class SemiMarkovData;
  class SemiMarkovIterator;
  class HiddenSemiMarkov;
  class NonhomogeneousMarkov;
  class NonhomogeneousMarkovData;
  class Sequences;
  class SequenceCharacteristics;

  class Switching_sequence;  // addition of Florence Chaubert

  /// \brief Categorical observation process for sequences

  class CategoricalSequenceProcess : public stat_tool::CategoricalProcess {  

  public :

    stat_tool::Distribution *length;  ///< sequence length distribution
    stat_tool::Curves *index_value;  ///< probabilities of each category as a function of the index parameter
    double *no_occurrence;  ///< probabilities of not observing each category
    stat_tool::Distribution **first_occurrence;  ///< time to the 1st occurrence distributions
    double *leave;          ///< probabilities of leaving definitively each category
    stat_tool::Distribution **recurrence_time;  ///< recurrence time distributions
    double *absorption;     ///< absorbing probabilities
    stat_tool::DiscreteParametric **sojourn_time;  ///< sojourn time distributions
    stat_tool::Distribution **nb_run;  ///< number of runs per sequence distributions
    stat_tool::Distribution **nb_occurrence;  ///< number of occurrences per sequence distributions

    void create_characteristic(const stat_tool::Distribution &ilength , bool* homogeneity ,
                               bool counting_flag = true);
    void create_characteristic(const stat_tool::Distribution &ilength , bool sojourn_time_flag = true ,
                               bool counting_flag = true);
    void copy(const CategoricalSequenceProcess &process , bool characteristic_flag = true);
    void init_occupancy(const CategoricalSequenceProcess &process , int occupancy_nb_value);
    void remove();

    CategoricalSequenceProcess(int inb_state = 0 , int inb_value = 0 ,
                               int observation_flag = false);
    CategoricalSequenceProcess(int inb_state , stat_tool::DiscreteParametric **occupancy);
    CategoricalSequenceProcess(const stat_tool::CategoricalProcess &process);
    CategoricalSequenceProcess(const CategoricalSequenceProcess &process ,
                               categorical_sequence_process_transformation transform = CATEGORICAL_SEQUENCE_PROCESS_COPY ,
                               int param = true);
    ~CategoricalSequenceProcess();
    CategoricalSequenceProcess& operator=(const CategoricalSequenceProcess &process);

    stat_tool::Distribution* weight_computation() const;

    static bool test_hidden(int nb_output_process , CategoricalSequenceProcess **process);

    static CategoricalSequenceProcess* occupancy_parsing(stat_tool::StatError &error , ifstream &in_file ,
                                                         int &line , const stat_tool::Chain &chain ,
                                                         double cumul_threshold = stat_tool::CUMUL_THRESHOLD);

    std::ostream& ascii_print(std::ostream &os , int process ,
                              stat_tool::FrequencyDistribution **empirical_observation ,
                              stat_tool::FrequencyDistribution *marginal_distribution ,
                              const SequenceCharacteristics *characteristics , bool exhaustive ,
                              bool file_flag , stat_tool::Forward **forward = NULL) const;
    std::ostream& spreadsheet_print(std::ostream &os , int process ,
                                    stat_tool::FrequencyDistribution **empirical_observation = NULL ,
                                    stat_tool::FrequencyDistribution *marginal_distribution = NULL ,
                                    const SequenceCharacteristics *characteristics = NULL ,
                                    stat_tool::Forward **forward = NULL) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    stat_tool::FrequencyDistribution **empirical_observation = NULL ,
                    stat_tool::FrequencyDistribution *marginal_distribution = NULL ,
                    const SequenceCharacteristics *characteristics = NULL ,
                    const stat_tool::FrequencyDistribution *length_distribution = NULL ,
                    stat_tool::Forward **forward = NULL) const;
    void plotable_write(stat_tool::MultiPlotSet &plot , int &index , int process ,
                        stat_tool::FrequencyDistribution **empirical_observation = NULL ,
                        stat_tool::FrequencyDistribution *marginal_distribution = NULL ,
                        const SequenceCharacteristics *characteristics = NULL ,
                        const stat_tool::FrequencyDistribution *length_distribution = NULL ,
                        stat_tool::Forward **forward = NULL) const;
  };


  /// \brief Correlation functions

  class Correlation : public stat_tool::StatInterface , public stat_tool::Curves {

    friend class VariableOrderMarkovChain;
    friend class VariableOrderMarkov;
    friend class Sequences;

    friend std::ostream& operator<<(std::ostream &os , const Correlation &correlation)
    { return correlation.ascii_write(os); }

  private :

    stat_tool::correlation_type type;  ///<  correlation coefficient type (PEARSON/SPEARMAN/KENDALL)
    correlation_variable_type *variable_type;  ///< variable types (OBSERVED/THEORETICAL STATE/OUTPUT)
    int *variable1;         ///< 1st variables
    int *variable2;         ///< 2nd variables
    autocorrelation_function_type function_type;  ///< theoretical correlation function type (AUTOREGRESSIVE/WHITE_NOISE)
    double *theoretical_function;  ///< theoretical correlation function for a first-order autoregressive model or a filtered white noise

    void copy(const Correlation &correl);
    void remove();
    bool plot_print(const char *path , double *confidence_limit) const;

  public :

    Correlation();
    Correlation(stat_tool::correlation_type itype , int max_lag , int ivariable1 , int ivariable2);
    Correlation(int inb_curve , int ilength , bool frequency_flag , stat_tool::correlation_type itype);
    Correlation(const Correlation &correl)
    :stat_tool::Curves(correl) { copy(correl); }
    ~Correlation();
    Correlation& operator=(const Correlation &correl);

    Correlation* merge(stat_tool::StatError &error , int nb_correl , const Correlation **icorrel) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    bool autoregressive_model_autocorrelation(stat_tool::StatError &error , double autoregressive_coeff);
    bool white_noise_correlation(stat_tool::StatError &error , int nb_point , double *filter ,
                                 int residual = true);
    bool white_noise_correlation(stat_tool::StatError &error , const stat_tool::Distribution &dist);
    bool white_noise_correlation(stat_tool::StatError &error , int order);

    // class member access

    stat_tool::correlation_type get_type() const { return type; }
    correlation_variable_type get_variable_type(int index) const { return variable_type[index]; }
    int get_variable1(int index) const { return variable1[index]; }
    int get_variable2(int index) const { return variable2[index]; }
    double get_theoretical_function(int lag) const { return theoretical_function[lag]; }
  };


  class MarkovianSequences;
  class TimeEvents;
  class RenewalData;

  /// \brief Sequences

  class Sequences : public stat_tool::StatInterface {

    friend std::ostream& operator<<(std::ostream &os , const Sequences &seq)
    { return seq.ascii_write(os); }

  protected :

    int nb_sequence;        ///< number of sequences
    int *identifier;        ///< sequence identifiers
    int max_length;         ///< maximum sequence length
    int cumul_length;       ///< cumulative sequence length
    int *length;            ///< sequence lengths
    stat_tool::FrequencyDistribution *length_distribution;  ///< sequence length frequency distribution
    int **vertex_identifier;  ///< identifiers of vertices of an associated MTG
    index_parameter_type index_param_type;  ///< index parameter type (TIME/POSITION)
    stat_tool::FrequencyDistribution *index_parameter_distribution;  ///< explicit index parameter frequency distribution
    stat_tool::FrequencyDistribution *index_interval;  ///< frequency distribution of the intervals between explicit index parameters
    int **index_parameter;  ///< explicit index parameters
    int nb_variable;        ///< number of variables
    stat_tool::variable_nature *type;  ///< variable types (INT_VALUE/REAL_VALUE/STATE)
    double *min_value;      ///<  minimum value for each variable
    double *max_value;      ///<  maximum value for each variable
    stat_tool::FrequencyDistribution **marginal_distribution;  ///< marginal frequency distributions
    stat_tool::Histogram **marginal_histogram;  ///< marginal histograms
    int ***int_sequence;    ///< sequences, integer-valued variables
    double ***real_sequence;  ///< sequences, real-valued variables

    void init(int inb_sequence , int *iidentifier , int *ilength , int **ivertex_identifier ,
              index_parameter_type iindex_param_type , int inb_variable , stat_tool::variable_nature *itype ,
              bool vertex_identifier_copy , bool init_flag);
    void init(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag);
    void copy(const Sequences &seq);
    void reverse(const Sequences &seq);
    void add_state_variable(const Sequences &seq);
    void explicit_index_parameter(const Sequences &seq);
    void remove_index_parameter(const Sequences &seq);
    void remove();

    bool increasing_index_parameter_checking(stat_tool::StatError &error , bool strict ,
                                             const char *pattern_label) const;
    bool increasing_sequence_checking(stat_tool::StatError &error , int variable , bool strict ,
                                      const char *pattern_label , const char *variable_label) const;

    void cluster(const Sequences &seq , int variable , int step , stat_tool::rounding mode);
    void transcode(const Sequences &seq , int ivariable , int min_category , int max_category ,
                   int *category , bool add_variable = false);
    void cluster(const Sequences &seq , int variable , int nb_class , double *limit);
    void select_variable(const Sequences &seq , int *variable);

    bool pointwise_average_ascii_print(stat_tool::StatError &error , const std::string path , int *frequency ,
                                       bool dispersion , sequence_type output) const;
    bool pointwise_average_spreadsheet_print(stat_tool::StatError &error , const std::string path , int *frequency ,
                                             bool dispersion , sequence_type output) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    std::ostream& ascii_print(std::ostream &os , output_sequence_format format , bool comment_flag ,
                              double *posterior_probability = NULL , double *entropy = NULL ,
                              double *nb_state_sequence = NULL , double *posterior_state_probability = NULL ,
                              int line_nb_character = stat_tool::LINE_NB_CHARACTER) const;
    bool plot_print(const char *path , int ilength) const;

    void max_length_computation();
    void cumul_length_computation();
    void build_length_frequency_distribution();

    void index_parameter_computation();
    void build_index_parameter_frequency_distribution();
    void index_interval_computation();

    void min_value_computation(int variable);
    void max_value_computation(int variable);
    void build_marginal_frequency_distribution(int variable);
    void build_marginal_histogram(int variable , double bin_width = stat_tool::D_DEFAULT ,
                                  double imin_value = stat_tool::D_INF);

    std::ostream& alignment_ascii_print(std::ostream &os , int width , int ref_index , int test_index ,
                                        const Sequences &alignment , int alignment_index) const;
    std::ostream& alignment_spreadsheet_print(std::ostream &os , int ref_index , int test_index ,
                                              const Sequences &alignment , int alignment_index) const;

    double indel_distance_computation(const stat_tool::VectorDistance &vector_dist ,
                                      double **rank , double **max_category_distance) const;
    double indel_distance_computation(const stat_tool::VectorDistance &vector_dist ,
                                      int index , int position , double **rank ,
                                      double **max_category_distance) const;
    double substitution_distance_computation(const stat_tool::VectorDistance &vector_dist , int ref_index ,
                                             int test_index , int ref_position , int test_position ,
                                             double **rank , const Sequences *test_seq = NULL) const;
    double substitution_distance_computation(int ref_index , int test_index , int ref_position,
                                             int test_position , double substitution_distance) const;

    std::ostream& multiple_alignment_ascii_print(std::ostream &os) const;
    bool multiple_alignment_ascii_print(stat_tool::StatError &error , const std::string path) const;

    Sequences* multiple_alignment(const Sequences &test_seq , const stat_tool::VectorDistance &vector_dist ,
                                  double **rank , double **max_category_distance , bool begin_free ,
                                  bool end_free , insertion_deletion_cost indel_cost , double indel_factor) const;

    void correlation_computation(Correlation &correl , int variable1 , int variable2 ,
                                 correlation_normalization normalization , bool individual_mean = false) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_segment ,
                                      double **profiles , const char *label ,
                                      double **piecewise_function = NULL , long double **change_point = NULL ,
                                      stat_tool::Distribution **segment_length = NULL ,
                                      stat_tool::Distribution *prior_segment_length = NULL ,
                                      long double **begin_conditonal_entropy = NULL ,
                                      long double **end_conditional_entropy = NULL ,
                                      long double **change_point_entropy = NULL) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_segment ,
                                            double **profiles , const char *label ,
                                            bool common_contrast = true , double ***piecewise_function = NULL ,
                                            long double **change_point = NULL ,
                                            stat_tool::Distribution **segment_length = NULL ,
                                            stat_tool::Distribution *prior_segment_length = NULL ,
                                            long double **begin_conditonal_entropy = NULL ,
                                            long double **end_conditional_entropy = NULL ,
                                            long double **change_point_entropy = NULL) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_segment ,
                                     double **profiles , bool common_contrast = true ,
                                     double ***piecewise_function = NULL ,
                                     long double **change_point = NULL ,
                                     stat_tool::Distribution **segment_length = NULL ,
                                     stat_tool::Distribution *prior_segment_length = NULL ,
                                     long double **begin_conditonal_entropy = NULL ,
                                     long double **end_conditional_entropy = NULL ,
                                     long double **change_point_entropy = NULL) const;
    void change_point_profile_plotable_write(stat_tool::MultiPlot &plot , int index , int nb_segment ,
                                             long double **change_point) const;
    void entropy_profile_plotable_write(stat_tool::MultiPlot &plot , int index ,
                                        long double *begin_conditional_entropy ,
                                        long double *end_conditional_entropy ,
                                        long double *change_point_entropy) const;

    void gamma_hyperparameter_computation(int index , int variable ,
                                          double *hyperparam) const;
    void gaussian_gamma_hyperparameter_computation(int index , int variable ,
                                                   double *hyperparam) const;
    int nb_parameter_computation(int index , int nb_segment , segment_model *model_type ,
                                 bool common_contrast) const;
    double one_segment_likelihood(int index , segment_model *model_type , bool common_contrast ,
                                  double *shape_parameter , double **rank) const;
    double piecewise_linear_function(int index , int variable , int nb_segment , segment_model model_type ,
                                     bool common_contrast , int *change_point , int *seq_index_parameter ,
                                     double **piecewise_function , double **imean = NULL , double **variance = NULL ,
                                     double *global_variance = NULL , double **iintercept = NULL , double **islope = NULL ,
                                     double **iautoregressive_coeff = NULL , double **correlation = NULL ,
                                     double **slope_standard_deviation = NULL , double **iindex_parameter_mean = NULL ,
                                     long double **iindex_parameter_variance = NULL , double **determination_coeff = NULL) const;
    std::ostream& piecewise_linear_function_ascii_print(std::ostream &os , int index , int variable , int nb_segment ,
                                                        segment_model model_type , bool common_contrast , int *change_point ,
                                                        int *seq_index_parameter , double **mean , double **variance ,
                                                        double **intercept , double **slope , double **autoregressive_coeff ,
                                                        double **correlation = NULL , double **slope_standard_deviation = NULL ,
                                                        double **index_parameter_mean = NULL , long double **index_parameter_variance = NULL ,
                                                        double **determination_coeff = NULL) const;
    std::ostream& piecewise_linear_function_spreadsheet_print(std::ostream &os , int index , int variable , int nb_segment ,
                                                              segment_model model_type , bool common_contrast , int *change_point ,
                                                              int *seq_index_parameter , double **mean , double **variance ,
                                                              double **intercept , double **slope , double **autoregressive_coeff ,
                                                              double **correlation = NULL , double **slope_standard_deviation = NULL ,
                                                              double **index_parameter_mean = NULL , long double **index_parameter_variance = NULL ,
                                                              double **determination_coeff = NULL) const;
    double continuous_piecewise_linear_function(std::ostream &os , int index , int variable , int nb_segment ,
                                                segment_model model_type , bool common_contrast ,
                                                int *change_point , int *seq_index_parameter ,
                                                double *intercept , double *slope ,
                                                double *corrected_intercept , double *corrected_slope) const;
    Sequences* segmentation_output(int nb_segment , segment_model *model_type , bool common_contrast , 
                                   bool display , sequence_type output = SEQUENCE , int *ichange_point = NULL ,
                                   bool continuity = false);
    void forward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                          double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                          double **seq_mean , int *seq_index_parameter , double **hyperparam ,
                          double **rank , long double *contrast , int nb_segment = 0) const;
    void backward_contrast(int time , int index , segment_model *model_type , bool common_contrast ,
                           double ***factorial , double *shape_parameter , double ***binomial_coeff ,
                           double **seq_mean , int *seq_index_parameter , double **hyperparam ,
                           double **rank , long double *contrast) const;
    double segmentation(int index , int nb_segment , segment_model *model_type ,
                        bool common_contrast , double *shape_parameter , double **rank ,
                        double *isegmentation_likelihood = NULL , int *nb_parameter = NULL ,
                        double *segment_penalty = NULL);
    double forward_backward(int index , int nb_segment , segment_model *model_type ,
                            bool common_contrast , double *shape_parameter , double **rank ,
                            double *likelihood , long double *segmentation_entropy ,
                            long double *first_order_entropy , long double *change_point_entropy ,
                            double *uniform_entropy , long double *marginal_entropy) const;
    double forward_backward(int index , int nb_segment , segment_model *model_type ,
                            bool common_contrast , double *shape_parameter ,  double **rank ,
                            std::ostream *os , stat_tool::MultiPlotSet *plot_set ,
                            double &segment_length_max , change_point_profile output ,
                            stat_tool::output_format format) const;
    double forward_backward_sampling(int index , int nb_segment , segment_model *model_type ,
                                     bool common_contrast , double *shape_parameter , double **rank ,
                                     std::ostream &os , stat_tool::output_format format ,
                                     int nb_segmentation) const;
    double N_segmentation(int index , int nb_segment , segment_model *model_type ,
                          bool common_contrast , double *shape_parameter , double **irank ,
                          std::ostream &os , stat_tool::output_format format , int inb_segmentation ,
                          double likelihood) const;
    double forward_backward_dynamic_programming(int index , int nb_segment , segment_model *model_type ,
                                                bool common_contrast , double *shape_parameter ,
                                                double **rank , std::ostream *os ,
                                                stat_tool::MultiPlotSet *plot_set , change_point_profile output ,
                                                stat_tool::output_format format , double likelihood = stat_tool::D_INF) const;

    std::ostream& profile_ascii_print(std::ostream &os , int index , int nb_state ,
                                      double **profiles , double *begin_conditional_entropy ,
                                      double *marginal_entropy , double *begin_partial_entropy ,
                                      double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    std::ostream& profile_spreadsheet_print(std::ostream &os , int index , int nb_state ,
                                            double **profiles , double *begin_conditional_entropy ,
                                            double *marginal_entropy , double *begin_partial_entropy ,
                                            double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    std::ostream& profile_plot_print(std::ostream &os , int index , int nb_state ,
                                     double **profiles , double *begin_conditional_entropy ,
                                     double *marginal_entropy , double *begin_partial_entropy ,
                                     double *end_conditional_entropy = NULL , double *end_partial_entropy = NULL) const;
    void profile_plotable_write(stat_tool::MultiPlot &plot , int index , int nb_state ,
                                double **profiles) const;
    void segment_length_distribution_plotable_write(stat_tool::MultiPlot &plot , int nb_segment ,
                                                    double segment_length_max ,
                                                    stat_tool::Distribution **segment_length ,
                                                    stat_tool::Distribution *prior_segment_length) const;
    void entropy_profile_plotable_write(stat_tool::MultiPlot &plot , int index , double *begin_entropy ,
                                        double *end_entropy = NULL , double *marginal_entropy = NULL) const;

    bool segment_profile_write(stat_tool::StatError &error , std::ostream &os , int iidentifier ,
                               int nb_segment , segment_model *model_type ,
                               bool common_contrast , double *shape_parameter ,
                               change_point_profile output = SEGMENT ,
                               stat_tool::output_format format = stat_tool::ASCII ,
                               stat_tool::latent_structure_algorithm segmentation = stat_tool::FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;

  public :

    Sequences();
    Sequences(int inb_sequence , int inb_variable);
    Sequences(int inb_sequence , int *iidentifier , int *ilength ,
              int **ivertex_identifier , index_parameter_type iindex_param_type , int inb_variable ,
              stat_tool::variable_nature *itype , bool vertex_identifier_copy = true , bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , ivertex_identifier ,
           iindex_param_type , inb_variable , itype ,
           vertex_identifier_copy , init_flag); }
    Sequences(int inb_sequence , int *iidentifier , int *ilength , index_parameter_type iindex_param_type ,  // AML interface
              int inb_variable , stat_tool::variable_nature itype , int ***iint_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,  // AML interface
              double ***ireal_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int **ivertex_identifier ,
              index_parameter_type iindex_param_type , int **iindex_parameter , int inb_variable ,
              stat_tool::variable_nature *itype , int ***iint_sequence , double ***ireal_sequence);
    Sequences(int inb_sequence , int *iidentifier , int *ilength , int inb_variable ,
              bool init_flag = false)
    { init(inb_sequence , iidentifier , ilength , inb_variable , init_flag); }
    Sequences(const stat_tool::FrequencyDistribution &ilength_distribution , int inb_variable ,
              stat_tool::variable_nature *itype , bool init_flag = false);
    Sequences(const RenewalData &timev);
    Sequences(const Sequences &seq , int variable , stat_tool::variable_nature itype);
    Sequences(const Sequences &seq , int inb_sequence , int *index);
    Sequences(const Sequences &seq , bool *auxiliary);
    Sequences(const Sequences &seq , sequence_transformation transform = SEQUENCE_COPY);
    ~Sequences();
    Sequences& operator=(const Sequences &seq);

    stat_tool::DiscreteDistributionData* extract(stat_tool::StatError &error , int variable) const;

    stat_tool::Vectors* build_vectors(bool index_variable) const;
    stat_tool::Vectors* extract_vectors(stat_tool::StatError &error , sequence_pattern pattern ,
                                        int variable = stat_tool::I_DEFAULT ,
                                        int value = stat_tool::I_DEFAULT) const;

    MarkovianSequences* markovian_sequences(stat_tool::StatError &error) const;

    bool check(stat_tool::StatError &error , const char *pattern_label);

    TimeEvents* extract_time_events(stat_tool::StatError &error , int variable ,
                                    int begin_date , int end_date , int previous_date = stat_tool::I_DEFAULT ,
                                    int next_date = stat_tool::I_DEFAULT) const;
    RenewalData* extract_renewal_data(stat_tool::StatError &error , int variable ,
                                      int begin_index_parameter , int end_index_parameter) const;

    Sequences* merge(stat_tool::StatError &error , int nb_sample , const Sequences **iseq) const;
    Sequences* merge(stat_tool::StatError &error , int nb_sample , const std::vector<Sequences> iseq) const;

    Sequences* shift(stat_tool::StatError &error , int variable , int shift_param) const;
    Sequences* shift(stat_tool::StatError &error , int variable , double shift_param) const;
    Sequences* thresholding(stat_tool::StatError &error , int variable , int threshold ,
                            stat_tool::threshold_direction mode) const;
    Sequences* thresholding(stat_tool::StatError &error , int variable , double threshold ,
                            stat_tool::threshold_direction mode) const;
    Sequences* cluster(stat_tool::StatError &error , int variable , int step ,
                       stat_tool::rounding mode = stat_tool::FLOOR) const;
    Sequences* transcode(stat_tool::StatError &error , int variable , int *category) const;
    Sequences* transcode(stat_tool::StatError &error , int variable , std::vector<int> category) const;
    Sequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                       int *ilimit) const;
    Sequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                       std::vector<int> ilimit) const;
    Sequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                       double *ilimit) const;
    Sequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                       std::vector<double> ilimit) const;
    Sequences* scaling(stat_tool::StatError &error , int variable , int scaling_coeff) const;
    Sequences* scaling(stat_tool::StatError &error , int variable , double scaling_coeff) const;
    Sequences* round(stat_tool::StatError &error , int variable = stat_tool::I_DEFAULT ,
                     stat_tool::rounding mode = stat_tool::ROUND) const;

    Sequences* index_parameter_select(stat_tool::StatError &error , bool display ,
                                      int min_index_parameter ,
                                      int max_index_parameter , bool keep) const;
    Sequences* value_select(stat_tool::StatError &error , bool display , int variable ,
                            int imin_value , int imax_value , bool keep = true) const;
    Sequences* value_select(stat_tool::StatError &error , bool display , int variable ,
                            double imin_value , double imax_value , bool keep = true) const;
    Sequences* select_individual(stat_tool::StatError &error , int inb_sequence , int *iidentifier ,
                                 bool keep = true) const;
    Sequences* select_individual(stat_tool::StatError &error , int inb_sequence , std::vector<int> iidentifier ,
                                 bool keep = true) const;

    Sequences* remove_index_parameter(stat_tool::StatError &error) const;
    Sequences* explicit_index_parameter(stat_tool::StatError &error) const;
    Sequences* select_variable(stat_tool::StatError &error , int inb_variable , int *ivariable ,
                               bool keep = true) const;
    Sequences* select_variable(stat_tool::StatError &error , int inb_variable , std::vector<int> ivariable ,
                               bool keep = true) const;
    Sequences* sum_variable(stat_tool::StatError &error , int nb_summed_variable , int *ivariable) const;
    Sequences* sum_variable(stat_tool::StatError &error , int nb_summed_variable , std::vector<int> ivariable) const;
    Sequences* merge_variable(stat_tool::StatError &error , int nb_sample , const Sequences **iseq ,
                              int ref_sample = stat_tool::I_DEFAULT) const;
    Sequences* merge_variable(stat_tool::StatError &error , int nb_sample , const std::vector<Sequences> iseq ,
                              int ref_sample = stat_tool::I_DEFAULT) const;
    Sequences* difference_variable(stat_tool::StatError &error , const Sequences &residual) const;
    Sequences* shift_variable(stat_tool::StatError &error , int variable , int lag) const;

    Sequences* reverse(stat_tool::StatError &error) const;
    Sequences* length_select(stat_tool::StatError &error , bool display , int min_length ,
                             int imax_length , bool keep = true) const;
    Sequences* remove_run(stat_tool::StatError &error , int variable , int ivalue ,
                          run_position position , int max_run_length = stat_tool::I_DEFAULT) const;
    Sequences* truncate(stat_tool::StatError &error , int max_index_parameter) const;
    Sequences* index_parameter_extract(stat_tool::StatError &error , int min_index_parameter ,
                                       int max_index_parameter = stat_tool::I_DEFAULT) const;
    Sequences* segmentation_extract(stat_tool::StatError &error , int variable , int nb_value ,
                                    int *ivalue , bool keep = true ,
                                    bool concatenation = false) const;
    Sequences* segmentation_extract(stat_tool::StatError &error , int variable , int nb_value ,
                                    std::vector<int> ivalue , bool keep = true ,
                                    bool concatenation = false) const;

    Sequences* cumulate(stat_tool::StatError &error , int variable = stat_tool::I_DEFAULT) const;
    Sequences* difference(stat_tool::StatError &error , int variable = stat_tool::I_DEFAULT ,
                          bool first_element = false) const;
    Sequences* log_transform(stat_tool::StatError &error , int variable = stat_tool::I_DEFAULT ,
                             stat_tool::log_base base = stat_tool::NATURAL) const;
    Sequences* relative_growth_rate(stat_tool::StatError &error , double growth_factor = GROWTH_FACTOR) const;
    Sequences* sequence_normalization(stat_tool::StatError &error , int variable = stat_tool::I_DEFAULT) const;
    Sequences* moving_average(stat_tool::StatError &error , int nb_point , double *filter ,
                              int variable = stat_tool::I_DEFAULT , bool begin_end = false ,
                              bool segmentation = false , sequence_type output = TREND) const;
    Sequences* moving_average(stat_tool::StatError &error , int nb_point , std::vector<double> filter ,
                              int variable = stat_tool::I_DEFAULT , bool begin_end = false ,
                              bool segmentation = false , sequence_type output = TREND) const;
    Sequences* moving_average(stat_tool::StatError &error , const stat_tool::Distribution &dist ,
                              int variable = stat_tool::I_DEFAULT , bool begin_end = false ,
                              bool segmentation = false , sequence_type output = TREND) const;

    Sequences* pointwise_average(stat_tool::StatError &error , bool robust = false , bool circular = false ,
                                 bool dispersion = false , sequence_type output = SEQUENCE ,
                                 const std::string path = "" ,
                                 stat_tool::output_format format = stat_tool::ASCII) const;

    bool mean_error_computation(stat_tool::StatError &error , bool display , int variable ,
                                int iidentifier = stat_tool::I_DEFAULT , bool robust = false) const;

    Sequences* recurrence_time_sequences(stat_tool::StatError &error , int variable , int value) const;
    Sequences* sojourn_time_sequences(stat_tool::StatError &error , int variable) const;

    Sequences* transform_position(stat_tool::StatError &error , int step) const;

    Sequences* cross(stat_tool::StatError &error) const;

    static Sequences* ascii_read(stat_tool::StatError &error , const std::string path ,
                                 bool old_format = false);

    std::ostream& line_write(std::ostream &os) const;

    virtual std::ostream& ascii_data_write(std::ostream &os , output_sequence_format format = COLUMN ,
                                           bool exhaustive = false) const;
    virtual bool ascii_data_write(stat_tool::StatError &error , const std::string path ,
                                  output_sequence_format format = COLUMN , bool exhaustive = false) const;
    bool plot_data_write(stat_tool::StatError &error , const char *prefix ,
                         const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable_data(stat_tool::StatError &error) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    int min_index_parameter_computation() const;
    int max_index_parameter_computation(bool last_position = false) const;

    void marginal_frequency_distribution_computation(int variable);
    bool select_bin_width(stat_tool::StatError &error , int variable , double bin_width ,
                          double imin_value = stat_tool::D_INF);

    double mean_computation(int variable) const;
    double variance_computation(int variable , double mean) const;
    double mean_absolute_deviation_computation(int variable , double location) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable , double mean , double variance) const;
    double kurtosis_computation(int variable , double mean , double variance) const;
    double* mean_direction_computation(int variable , stat_tool::angle_unit unit) const;

    stat_tool::FrequencyDistribution* value_index_interval_computation(stat_tool::StatError &error ,
                                                                       int variable , int value) const;

    Correlation* correlation_computation(stat_tool::StatError &error , int variable1 , int variable2 ,
                                         stat_tool::correlation_type itype = stat_tool::PEARSON ,
                                         int max_lag = stat_tool::I_DEFAULT ,
                                         correlation_normalization normalization = EXACT ,
                                         bool individual_mean = false) const;
    Correlation* partial_autocorrelation_computation(stat_tool::StatError &error , int variable ,
                                                     stat_tool::correlation_type itype = stat_tool::PEARSON ,
                                                     int max_lag = stat_tool::I_DEFAULT) const;

    stat_tool::DistanceMatrix* alignment(stat_tool::StatError &error , bool display ,
                                         const stat_tool::VectorDistance &ivector_dist ,
                                         int ref_identifier = stat_tool::I_DEFAULT , int test_identifier = stat_tool::I_DEFAULT ,
                                         bool begin_free = false , bool end_free = false ,
                                         insertion_deletion_cost indel_cost = ADAPTATIVE ,
                                         double indel_factor = INDEL_FACTOR_1 , bool transposition_flag = false ,
                                         double transposition_factor = TRANSPOSITION_FACTOR ,
                                         const std::string result_path = "" , stat_tool::output_format result_format = stat_tool::ASCII ,
                                         const std::string alignment_path = "") const;
    stat_tool::DistanceMatrix* alignment(stat_tool::StatError &error , bool display ,
                                         int ref_identifier = stat_tool::I_DEFAULT , int test_identifier = stat_tool::I_DEFAULT ,
                                         bool begin_free = false , bool end_free = false ,
                                         const std::string result_path = "" , stat_tool::output_format result_format = stat_tool::ASCII ,
                                         const std::string alignment_path = "") const;

    Sequences* multiple_alignment(stat_tool::StatError &error , bool display ,
                                  const stat_tool::VectorDistance &ivector_dist ,
                                  bool begin_free = false , bool end_free = false ,
                                  insertion_deletion_cost indel_cost = ADAPTATIVE ,
                                  double indel_factor = INDEL_FACTOR_N ,
                                  stat_tool::hierarchical_strategy strategy = stat_tool::AGGLOMERATIVE ,
                                  const std::string path = "") const;

    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int nb_segment , int *ichange_point , segment_model *model_type ,
                            bool common_contrast , double *shape_parameter ,
                            sequence_type output = SEQUENCE , bool continuity = false) const;
    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int nb_segment , std::vector<int> ichange_point , std::vector<segment_model> model_type ,
                            bool common_contrast , std::vector<double> shape_parameter ,
                            sequence_type output = SEQUENCE , bool continuity = false) const;
    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int nb_segment , segment_model *model_type ,
                            bool common_contrast , double *shape_parameter ,
                            sequence_type output , bool continuity = false) const;
    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int nb_segment , std::vector<segment_model> model_type ,
                            bool common_contrast , std::vector<double> shape_parameter ,
                            sequence_type output , bool continuity = false) const;
    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int max_nb_segment , segment_model *model_type ,
                            bool common_contrast , double *shape_parameter ,
                            stat_tool::model_selection_criterion criterion = stat_tool::LIKELIHOOD_SLOPE ,
                            int min_nb_segment = 0 , int penalty_shape_type = 2 ,
                            sequence_type output = SEQUENCE) const;
    Sequences* segmentation(stat_tool::StatError &error , bool display , int iidentifier ,
                            int max_nb_segment , std::vector<segment_model> model_type ,
                            bool common_contrast , std::vector<double> shape_parameter ,
                            stat_tool::model_selection_criterion criterion = stat_tool::LIKELIHOOD_SLOPE ,
                            int min_nb_segment = 0 , int penalty_shape_type = 2 ,
                            sequence_type output = SEQUENCE) const;

//    Sequences* hierarchical_segmentation(stat_tool::StatError &error , std::ostream &os , int iidentifier ,
//                                         int max_nb_segment , segment_model *model_type) const;

    bool segment_profile_ascii_write(stat_tool::StatError &error , int iidentifier ,
                                     int nb_segment , std::vector<segment_model> model_type ,
                                     bool common_contrast , std::vector<double> shape_parameter ,
                                     change_point_profile output = SEGMENT ,
                                     stat_tool::latent_structure_algorithm segmentation = stat_tool::FORWARD_DYNAMIC_PROGRAMMING ,
                                     int nb_segmentation = NB_SEGMENTATION) const;
    bool segment_profile_write(stat_tool::StatError &error , const std::string path , int iidentifier ,
                               int nb_segment , std::vector<segment_model> model_type ,
                               bool common_contrast , std::vector<double> shape_parameter ,
                               change_point_profile output = SEGMENT ,
                               stat_tool::output_format format = stat_tool::ASCII ,
                               stat_tool::latent_structure_algorithm segmentation = stat_tool::FORWARD_DYNAMIC_PROGRAMMING ,
                               int nb_segmentation = NB_SEGMENTATION) const;

    bool segment_profile_plot_write(stat_tool::StatError &error , const char *prefix ,
                                    int iidentifier , int nb_segment , segment_model *model_type ,
                                    bool common_contrast , double *shape_parameter ,
                                    change_point_profile output = SEGMENT , const char *title = NULL) const;
    stat_tool::MultiPlotSet* segment_profile_plotable_write(stat_tool::StatError &error , int iidentifier ,
                                                            int nb_segment , segment_model *model_type ,
                                                            bool common_contrast , double *shape_parameter ,
                                                            change_point_profile output = SEGMENT) const;
    stat_tool::MultiPlotSet* segment_profile_plotable_write(stat_tool::StatError &error , int iidentifier ,
                                                            int nb_segment , std::vector<segment_model> model_type ,
                                                            bool common_contrast , std::vector<double> shape_parameter ,
                                                            change_point_profile output = SEGMENT) const;

    // class member access

    int get_nb_sequence() const { return nb_sequence; }
    int get_identifier(int iseq) const { return identifier[iseq]; }
    int get_max_length() const { return max_length; }
    int get_cumul_length() const { return cumul_length; }
    int get_length(int index_seq) const { return length[index_seq]; }
    stat_tool::FrequencyDistribution* get_length_distribution() const { return length_distribution; }
    int get_vertex_identifier(int iseq , int index) const
    { return vertex_identifier[iseq][index]; }
    index_parameter_type get_index_param_type() const { return index_param_type; }
    stat_tool::FrequencyDistribution* get_index_parameter_distribution() const { return index_parameter_distribution; }
    stat_tool::FrequencyDistribution* get_index_interval() const { return index_interval; }
    int get_index_parameter(int iseq , int index) const
    { return index_parameter[iseq][index]; }
    int get_nb_variable() const { return nb_variable; }
    stat_tool::variable_nature get_type(int variable) const { return type[variable]; }
    double get_min_value(int variable) const { return min_value[variable]; }
    double get_max_value(int variable) const { return max_value[variable]; }
    stat_tool::FrequencyDistribution* get_marginal_distribution(int variable) const
    { return marginal_distribution[variable]; }
    stat_tool::Histogram* get_marginal_histogram(int variable) const
    { return marginal_histogram[variable]; }
    int get_int_sequence(int iseq , int variable , int index) const
    { return int_sequence[iseq][variable][index]; }
    double get_real_sequence(int iseq , int variable , int index) const
    { return real_sequence[iseq][variable][index]; }
  };


  /// \brief Sequence characteristics for a categorical variable

  class SequenceCharacteristics {

  public :

    int nb_value;           ///< number of categories
    stat_tool::Curves *index_value;    ///< empirical probabilities of each category as a function of the index parameter
    stat_tool::Curves *explicit_index_value;    ///< empirical probabilities of each category as a function of the explicit index parameter
    stat_tool::FrequencyDistribution **first_occurrence;  ///< time to the 1st occurrence frequency distributions
    stat_tool::FrequencyDistribution **recurrence_time;  ///< recurrence time frequency distributions
    stat_tool::FrequencyDistribution **sojourn_time;  ///< complete sojourn time frequency distributions
    stat_tool::FrequencyDistribution **initial_run;  ///< left-censored sojourn time frequency distributions
    stat_tool::FrequencyDistribution **final_run;  ///< right-censored sojourn time frequency distributions
    stat_tool::FrequencyDistribution **nb_run;  ///< number of runs per sequence frequency distributions
    stat_tool::FrequencyDistribution **nb_occurrence;  ///< number of occurrences per sequence frequency distributions

    void copy(const SequenceCharacteristics &characteristics);
    void reverse(const SequenceCharacteristics &characteristics);
    void remove();

    void create_sojourn_time_frequency_distribution(int max_length , int initial_run_flag = false);

    std::ostream& ascii_print(std::ostream &os , int type ,
                              const stat_tool::FrequencyDistribution &length_distribution ,
                              bool exhaustive , bool comment_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os , int type ,
                                    const stat_tool::FrequencyDistribution &length_distribution) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable , int type , const stat_tool::FrequencyDistribution &length_distribution) const;
    void plotable_write(stat_tool::MultiPlotSet &plot , int &index , int variable ,
                        int type , const stat_tool::FrequencyDistribution &length_distribution) const;

    SequenceCharacteristics(int inb_value = stat_tool::I_DEFAULT);
    SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                            bool initial_run_flag);
    SequenceCharacteristics(const SequenceCharacteristics &characteristics ,
                            sequence_transformation transform = SEQUENCE_COPY);
    ~SequenceCharacteristics();
    SequenceCharacteristics& operator=(const SequenceCharacteristics &characteristics);
  };


  class Function;

  /// \brief Self-transition probabilitiy as a function of the index parameter

  class SelfTransition : public stat_tool::Curves {

  public :

    SelfTransition(int ilength)
    :stat_tool::Curves(1 , ilength , true , false) {}

    Function* monomolecular_regression() const;
    Function* logistic_regression() const;
  };


  class VariableOrderMarkovChain;
  class VariableOrderMarkovChainData;

  /// \brief Sequences potentially generated by a (hidden) Markovian process

  class MarkovianSequences : public Sequences {

    friend class VariableOrderMarkovChain;
    friend class VariableOrderMarkov;
    friend class HiddenVariableOrderMarkov;
    friend class SemiMarkov;
    friend class HiddenSemiMarkov;
    friend class NonhomogeneousMarkov;

    friend std::ostream& operator<<(std::ostream &os , const MarkovianSequences &seq)
    { return seq.ascii_write(os); }

  protected :

    double *min_interval;   ///< minimum intervals between 2 values
    SelfTransition **self_transition;  ///< self transition probability as a function of the index parameter
    stat_tool::FrequencyDistribution ***observation_distribution;  ///< observation frequency distributions
    stat_tool::Histogram ***observation_histogram;  ///< observation histograms
    SequenceCharacteristics **characteristics;  ///< characteristics for categorical variables

    void init();
    void copy(const MarkovianSequences &seq , initial_run param = UNCHANGED);
    void reverse(const MarkovianSequences &seq);
    void add_state_variable(const MarkovianSequences &seq , initial_run param);
    void remove();

    MarkovianSequences* transcode(stat_tool::StatError &error ,
                                  const CategoricalSequenceProcess *process) const;
    MarkovianSequences* build_auxiliary_variable(stat_tool::DiscreteParametricProcess **discrete_process ,
                                                 stat_tool::ContinuousParametricProcess **continuous_process) const;
    MarkovianSequences* residual_sequences(CategoricalSequenceProcess **categorical_process ,
                                           stat_tool::DiscreteParametricProcess **discrete_process ,
                                           stat_tool::ContinuousParametricProcess **continuous_process) const;

    MarkovianSequences* remove_variable_1() const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    bool plot_print(const char *prefix , const char *title , int variable ,
                    int nb_variable) const;
    void plotable_write(stat_tool::MultiPlotSet &plot , int &index , int variable) const;

    void state_variable_init(stat_tool::variable_nature itype = stat_tool::STATE);

    void min_interval_computation(int variable);

    void transition_count_computation(const VariableOrderMarkovChainData &chain_data ,
                                      const VariableOrderMarkovChain &markov ,
                                      bool begin = true , bool non_terminal = false) const;
    void transition_count_computation(const stat_tool::ChainData &chain_data ,
                                      const SemiMarkov *smarkov = NULL) const;

    void self_transition_computation(int state);
    stat_tool::Distribution* weight_computation() const;
    void observation_frequency_distribution_computation(int variable , int nb_state);
    bool test_hidden(int variable) const;

    void build_index_value(int variable);
    void build_explicit_index_value(int variable);
    void build_first_occurrence_frequency_distribution(int variable);
    void build_recurrence_time_frequency_distribution(int variable);
    void build_sojourn_time_frequency_distribution(int variable , int initial_run_flag = false);
    void build_nb_run_frequency_distribution(int variable);
    void build_nb_occurrence_frequency_distribution(int variable);

    void censored_sojourn_time_frequency_distribution_computation(stat_tool::FrequencyDistribution **initial_run ,
                                                                  stat_tool::FrequencyDistribution **final_run ,
                                                                  stat_tool::FrequencyDistribution **single_run) const;

    std::ostream& linear_model_spreadsheet_print(std::ostream &os , int variable ,
                                                 stat_tool::ContinuousParametricProcess *process) const;
    bool linear_model_plot_print(const char *prefix , const char *title , int variable ,
                                 stat_tool::ContinuousParametricProcess *process) const;
    void linear_model_plotable_write(stat_tool::MultiPlotSet &plot , int &index , int variable ,
                                     stat_tool::ContinuousParametricProcess *process) const;

    void autocorrelation_computation(Correlation &correl , int state , int variable) const;
    std::ostream& autoregressive_model_ascii_print(std::ostream &os , int variable ,
                                                   stat_tool::ContinuousParametricProcess *process ,
                                                   bool file_flag) const;
    std::ostream& autoregressive_model_spreadsheet_print(std::ostream &os , int variable ,
                                                         stat_tool::ContinuousParametricProcess *process) const;
    bool autoregressive_model_plot_print(const char *prefix , const char *title , int variable ,
                                         stat_tool::ContinuousParametricProcess *process) const;
    void autoregressive_model_plotable_write(stat_tool::MultiPlotSet &plot , int &index , int variable ,
                                             stat_tool::ContinuousParametricProcess *process) const;

    template <typename Type>
    void gamma_estimation(Type ***state_sequence_count , int variable ,
                          stat_tool::ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void zero_inflated_gamma_estimation(Type ***state_sequence_count , int variable ,
                                        stat_tool::ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void inverse_gaussian_estimation(Type ***state_sequence_count , int variable ,
                                     stat_tool::ContinuousParametricProcess *process) const;
    template <typename Type>
    void gaussian_estimation(Type ***state_sequence_count , int variable ,
                             stat_tool::ContinuousParametricProcess *process) const;
    template <typename Type>
    void von_mises_estimation(Type ***state_sequence_count , int variable ,
                              stat_tool::ContinuousParametricProcess *process) const;
    template <typename Type>
    void  linear_model_estimation(Type ***state_sequence_count , int variable ,
                                  stat_tool::ContinuousParametricProcess *process) const;
    template <typename Type>
    void autoregressive_model_estimation(Type ***state_sequence_count , int variable ,
                                         stat_tool::ContinuousParametricProcess *process) const;

    std::ostream& likelihood_write(std::ostream &os , int nb_model , double **likelihood ,
                                   const char *label , bool exhaustive = false ,
                                   stat_tool::latent_structure_algorithm algorithm = stat_tool::NO_LATENT_STRUCTURE) const;
    bool likelihood_write(stat_tool::StatError &error , const std::string path , int nb_model ,
                          double **likelihood , const char *label ,
                          stat_tool::latent_structure_algorithm algorithm = stat_tool::NO_LATENT_STRUCTURE) const;

  public :

    MarkovianSequences();
    MarkovianSequences(int inb_sequence , int *iidentifier , int *ilength ,
                       int **ivertex_identifier , index_parameter_type iindex_param_type , int inb_variable ,
                       stat_tool::variable_nature *itype , bool vertex_identifier_copy = true , bool init_flag = false)
    :Sequences(inb_sequence , iidentifier , ilength , ivertex_identifier ,
               iindex_param_type , inb_variable , itype ,
               vertex_identifier_copy , init_flag) { init(); }
    MarkovianSequences(const stat_tool::FrequencyDistribution &ilength_distribution , int inb_variable ,
                       stat_tool::variable_nature *itype , bool init_flag = false)
    :Sequences(ilength_distribution , inb_variable , itype , init_flag) { init(); }
    MarkovianSequences(const MarkovianSequences &seq , int variable , stat_tool::variable_nature itype)
    :Sequences(seq , variable , itype) { init(); }
    MarkovianSequences(const Sequences &seq);
    MarkovianSequences(const MarkovianSequences &seq , bool *auxiliary);
    MarkovianSequences(const MarkovianSequences &seq , sequence_transformation transform = SEQUENCE_COPY ,
                       initial_run param = UNCHANGED);
    ~MarkovianSequences();
    MarkovianSequences& operator=(const MarkovianSequences &seq);

    stat_tool::DiscreteDistributionData* extract(stat_tool::StatError &error , stat_tool::process_distribution type ,
                                                 int variable , int value) const;

    MarkovianSequences* merge(stat_tool::StatError &error , int nb_sample ,
                              const MarkovianSequences **iseq) const;
    MarkovianSequences* merge(stat_tool::StatError &error , int nb_sample ,
                              const std::vector<MarkovianSequences> iseq) const;

    MarkovianSequences* cluster(stat_tool::StatError &error , int variable , int step ,
                                stat_tool::rounding mode = stat_tool::FLOOR) const;
    MarkovianSequences* transcode(stat_tool::StatError &error , int ivariable , int *category ,
                                  bool add_variable = false) const;
    MarkovianSequences* transcode(stat_tool::StatError &error , int ivariable , std::vector<int> category ,
                                  bool add_variable = false) const;
    MarkovianSequences* consecutive_values(stat_tool::StatError &error , bool display ,
                                           int ivariable , bool add_variable = false) const;
    MarkovianSequences* cluster(stat_tool::StatError &error , int ivariable , int nb_class ,
                                int *ilimit , bool add_variable = false) const;
    MarkovianSequences* cluster(stat_tool::StatError &error , int ivariable , int nb_class ,
                                std::vector<int> ilimit , bool add_variable = false) const;
    MarkovianSequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                                double *ilimit) const;
    MarkovianSequences* cluster(stat_tool::StatError &error , int variable , int nb_class ,
                                std::vector<double> ilimit) const;

    MarkovianSequences* remove_index_parameter(stat_tool::StatError &error) const;
    MarkovianSequences* explicit_index_parameter(stat_tool::StatError &error) const;
    MarkovianSequences* select_variable(stat_tool::StatError &error , int inb_variable ,
                                        int *ivariable , bool keep = true) const;
    MarkovianSequences* select_variable(stat_tool::StatError &error , int inb_variable ,
                                        std::vector<int> ivariable , bool keep = true) const;
    MarkovianSequences* merge_variable(stat_tool::StatError &error , int nb_sample ,
                                       const MarkovianSequences **iseq ,
                                       int ref_sample = stat_tool::I_DEFAULT) const;
    MarkovianSequences* merge_variable(stat_tool::StatError &error , int nb_sample ,
                                       const std::vector<MarkovianSequences> iseq ,
                                       int ref_sample = stat_tool::I_DEFAULT) const;

    MarkovianSequences* initial_run_computation(stat_tool::StatError &error) const;
    MarkovianSequences* add_absorbing_run(stat_tool::StatError &error ,
                                          int run_length = stat_tool::I_DEFAULT ,
                                          int sequence_length = stat_tool::I_DEFAULT ,
                                          bool add_variable = false) const;

    MarkovianSequences* split(stat_tool::StatError &error , int step) const;

    std::ostream& ascii_data_write(std::ostream &os , output_sequence_format format = COLUMN ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(stat_tool::StatError &error , const std::string path ,
                          output_sequence_format format = COLUMN , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const std::string path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    bool transition_count(stat_tool::StatError &error , bool display , int max_order ,
                          bool begin = false , transition_estimator estimator = MAXIMUM_LIKELIHOOD ,
                          const std::string path = "") const;
    bool word_count(stat_tool::StatError &error , bool display , int variable , int word_length ,
                    int begin_state = stat_tool::I_DEFAULT , int end_state = stat_tool::I_DEFAULT ,
                    int min_frequency = 1) const;
    bool mtg_write(stat_tool::StatError &error , const std::string path , stat_tool::variable_type *itype) const;

    int cumulative_distribution_function_computation(int variable , double **cdf) const;
    int cumulative_distribution_function_computation(int variable , int state , double **cdf) const;

    double iid_information_computation() const;

    void self_transition_computation();
    void self_transition_computation(bool *homogeneity);
    void sojourn_time_frequency_distribution_computation(int variable);

    void build_observation_frequency_distribution(int nb_state);
    void build_observation_histogram(int variable , int nb_state , double bin_width = stat_tool::D_DEFAULT);
    void build_observation_histogram(int nb_state);
    bool select_bin_width(stat_tool::StatError &error , int variable , double bin_width ,
                          double imin_value = stat_tool::D_INF);

    void build_characteristic(int variable = stat_tool::I_DEFAULT , bool sojourn_time_flag = true ,
                              bool initial_run_flag = false);

    NonhomogeneousMarkov* nonhomogeneous_markov_estimation(stat_tool::StatError &error , stat_tool::parametric_function *ident ,
                                                           bool counting_flag = true) const;

    VariableOrderMarkov* variable_order_markov_estimation(stat_tool::StatError &error , bool display ,
                                                          stat_tool::process_type itype , int min_order = 0 ,
                                                          int max_order = stat_tool::ORDER ,
                                                          memory_tree_selection algorithm = LOCAL_BIC ,
                                                          double threshold = LOCAL_BIC_THRESHOLD ,
                                                          transition_estimator estimator = LAPLACE ,
                                                          bool global_initial_transition = true ,
                                                          bool global_sample = true ,
                                                          bool counting_flag = true) const;
    VariableOrderMarkov* variable_order_markov_estimation(stat_tool::StatError &error ,
                                                          const VariableOrderMarkov &imarkov ,
                                                          bool global_initial_transition = true ,
                                                          bool counting_flag = true) const;
    VariableOrderMarkov* variable_order_markov_estimation(stat_tool::StatError &error ,
                                                          stat_tool::process_type itype , int order = 1 ,
                                                          bool global_initial_transition = true ,
                                                          bool counting_flag = true) const;

    VariableOrderMarkov* lumpability_estimation(stat_tool::StatError &error , bool display , int *category ,
                                                stat_tool::model_selection_criterion criterion = stat_tool::BIC ,
                                                int order = 1 , bool counting_flag = true) const;

    SemiMarkov* semi_markov_estimation(stat_tool::StatError &error , bool display , stat_tool::process_type itype ,
                                       stat_tool::censoring_estimator estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                                       bool counting_flag = true , int nb_iter = stat_tool::I_DEFAULT ,
                                       stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED) const;

    HiddenVariableOrderMarkov* hidden_variable_order_markov_estimation(stat_tool::StatError &error , bool display ,
                                                                       const HiddenVariableOrderMarkov &ihmarkov ,
                                                                       bool global_initial_transition = true ,
                                                                       bool common_dispersion = false ,
                                                                       bool counting_flag = true ,
                                                                       bool state_sequence = true ,
                                                                       int nb_iter = stat_tool::I_DEFAULT) const;
    HiddenVariableOrderMarkov* hidden_variable_order_markov_stochastic_estimation(stat_tool::StatError &error , bool display ,
                                                                                  const HiddenVariableOrderMarkov &ihmarkov ,
                                                                                  bool global_initial_transition = true ,
                                                                                  bool common_dispersion = false ,
                                                                                  int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                                  int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                                  double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                                  bool counting_flag = true ,
                                                                                  bool state_sequence = true ,
                                                                                  int nb_iter = stat_tool::I_DEFAULT) const;

    HiddenSemiMarkov* hidden_semi_markov_estimation(stat_tool::StatError &error , bool display ,
                                                    const HiddenSemiMarkov &ihsmarkov ,
                                                    bool poisson_geometric = false ,
                                                    bool common_dispersion = false ,
                                                    stat_tool::censoring_estimator estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                                                    bool counting_flag = true ,
                                                    bool state_sequence = true ,
                                                    int nb_iter = stat_tool::I_DEFAULT ,
                                                    stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED) const;
    HiddenSemiMarkov* hidden_semi_markov_estimation(stat_tool::StatError &error , bool display ,
                                                    stat_tool::process_type itype , int nb_state , bool left_right ,
                                                    double occupancy_mean = stat_tool::D_DEFAULT ,
                                                    bool poisson_geometric = false ,
                                                    bool common_dispersion = false ,
                                                    stat_tool::censoring_estimator estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                                                    bool counting_flag = true ,
                                                    bool state_sequence = true ,
                                                    int nb_iter = stat_tool::I_DEFAULT ,
                                                    stat_tool::duration_distribution_mean_estimator mean_estimator = stat_tool::COMPUTED) const;
    HiddenSemiMarkov* hidden_semi_markov_stochastic_estimation(stat_tool::StatError &error , bool display ,
                                                               const HiddenSemiMarkov &ihsmarkov ,
                                                               bool poisson_geometric = false ,
                                                               bool common_dispersion = false ,
                                                               int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                               int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                               double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                               stat_tool::censoring_estimator estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                                                               bool counting_flag = true ,
                                                               bool state_sequence = true ,
                                                               int nb_iter = stat_tool::I_DEFAULT) const;
    HiddenSemiMarkov* hidden_semi_markov_stochastic_estimation(stat_tool::StatError &error , bool display ,
                                                               stat_tool::process_type itype , int nb_state , bool left_right ,
                                                               double occupancy_mean = stat_tool::D_DEFAULT ,
                                                               bool poisson_geometric = false ,
                                                               bool common_dispersion = false ,
                                                               int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                               int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                               double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                               stat_tool::censoring_estimator estimator = stat_tool::COMPLETE_LIKELIHOOD ,
                                                               bool counting_flag = true ,
                                                               bool state_sequence = true ,
                                                               int nb_iter = stat_tool::I_DEFAULT) const;

    bool lumpability_test(stat_tool::StatError &error , bool display , int *category , int order = 1) const;

    bool comparison(stat_tool::StatError &error , bool display , int nb_model ,
                    const VariableOrderMarkov **imarkov , const std::string path = "") const;

    bool comparison(stat_tool::StatError &error , bool display , int nb_model ,
                    const SemiMarkov **ismarkov , const std::string path = "") const;

    bool comparison(stat_tool::StatError &error , bool display , int nb_model ,
                    const HiddenVariableOrderMarkov **ihmarkov ,
                    stat_tool::latent_structure_algorithm algorithm = stat_tool::FORWARD ,
                    const std::string path = "") const;

    bool comparison(stat_tool::StatError &error , bool display , int nb_model ,
                    const HiddenSemiMarkov **ihsmarkov ,
                    stat_tool::latent_structure_algorithm algorithm = stat_tool::FORWARD ,
                    const std::string path = "") const;

    // class member access

    stat_tool::Curves* get_self_transition(int state) const { return self_transition[state]; }
    stat_tool::FrequencyDistribution*** get_observation_distribution() const
    { return observation_distribution; }
    stat_tool::FrequencyDistribution** get_observation_distribution(int variable) const
    { return observation_distribution[variable]; }
    stat_tool::FrequencyDistribution* get_observation_distribution(int variable , int state) const
    { return observation_distribution[variable][state]; }
    stat_tool::Histogram*** get_observation_histogram() const { return observation_histogram; }
    stat_tool::Histogram** get_observation_histogram(int variable) const
    { return observation_histogram[variable]; }
    stat_tool::Histogram* get_observation_histogram(int variable , int state) const
    { return observation_histogram[variable][state]; }
    SequenceCharacteristics* get_characteristics(int variable) const
    { return characteristics[variable]; }
  };


};  // namespace sequence_analysis



#include "continuous_parametric_sequence_estimation.hpp"



#endif
