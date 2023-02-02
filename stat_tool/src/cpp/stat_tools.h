/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: stat_tools.h 18656 2015-11-03 07:25:48Z guedon $
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



#ifndef STAT_TOOLS_H
#define STAT_TOOLS_H


#include <fstream>
#include "reestimation.h"
#include "plotable.h"


namespace stat_tool {



/****************************************************************
 *
 *  Macros
 */


  #ifndef MIN
  #define MIN(x,y)  ((x) < (y) ? (x) : (y))
  #endif

  #ifndef MAX
  #define MAX(x,y)  ((x) > (y) ? (x) : (y))
  #endif



/****************************************************************
 *
 *  Constants
 */


  const int ERROR_LENGTH = 200;

  enum output_format {
    ASCII ,
    SPREADSHEET ,
//    BINARY ,
    GNUPLOT ,
    PLOT
  };

  const int I_DEFAULT = -1;              // default int
  const double D_DEFAULT = -1.;          // default double
  const double D_INF = -1.e37;           // smallest  real number
  const double DOUBLE_ERROR = 1.e-6;     // error on a sum of doubles
//  const double DOUBLE_ERROR = 5.e-6;      error on a sum of doubles

  enum test_distribution {
    STANDARD_NORMAL ,
    CHI2 ,
    FISHER ,
    STUDENT
  };

  const int NB_CRITICAL_PROBABILITY = 2;
  const double ref_critical_probability[NB_CRITICAL_PROBABILITY] = {0.05 , 0.01};

  const int NB_VALUE = 1000;             // number of values of a discrete variable
  const int SAMPLE_NB_VALUE = NB_VALUE;  // number of values of a discrete sample

  enum frequency_distribution_transformation {
    FREQUENCY_DISTRIBUTION_COPY ,
    SHIFT ,
    CLUSTER
  };

  enum rounding {
    FLOOR ,
    ROUND ,
    CEIL
  };

  enum discrete_parametric {
    CATEGORICAL ,
    BINOMIAL ,
    POISSON ,
    NEGATIVE_BINOMIAL ,
    POISSON_GEOMETRIC ,
    UNIFORM ,
    PRIOR_SEGMENT_LENGTH ,
    MULTINOMIAL                          // addition by Florence Chaubert
  };

  enum distribution_transformation {
    DISTRIBUTION_COPY ,
    NORMALIZATION
  };

  enum distribution_computation {
    STANDARD ,
    RENEWAL
  };

  enum continuous_parametric {
    GAMMA ,
    INVERSE_GAUSSIAN ,
    GAUSSIAN ,
    VON_MISES ,
    ZERO_INFLATED_GAMMA ,
    LINEAR_MODEL ,
    AUTOREGRESSIVE_MODEL
  };

  enum angle_unit {
    DEGREE ,
    RADIAN
  };

  enum compound_distribution {
    SUM ,
    ELEMENTARY ,
    COMPOUND
  };

  enum variable_type {
    NOMINAL ,
    ORDINAL ,
    NUMERIC ,
    CIRCULAR
  };

  enum log_base {
    NATURAL ,
    TWO ,
    TEN
  };

  enum process_type {
    ORDINARY ,
    EQUILIBRIUM ,
    DEFAULT_TYPE
  };

  enum estimation_criterion {
    LIKELIHOOD ,
    PENALIZED_LIKELIHOOD ,
    PARAMETRIC_REGULARIZATION
  };

  enum duration_distribution_mean_estimator {
    COMPUTED ,
    ESTIMATED ,
    ONE_STEP_LATE
  };

  enum censoring_estimator {
    PARTIAL_LIKELIHOOD ,
    COMPLETE_LIKELIHOOD ,
    KAPLAN_MEIER
  };

  enum model_selection_criterion {
    AIC ,
    AICc ,
    BIC ,
    BICc ,
    ICL ,
    ICLc ,
    mBIC ,
    LIKELIHOOD_SLOPE ,
    DIMENSION_JUMP ,
    SEGMENTATION_LIKELIHOOD_SLOPE ,
    SEGMENTATION_DIMENSION_JUMP
  };

  enum error_type {
    ERROR ,
    WARNING
  };

  enum variable_nature {
    INT_VALUE ,                          // integer-valued variable
    REAL_VALUE ,                         // real-valued variable
    STATE ,                              // state variable
    OLD_INT_VALUE ,                      // for ascending compatibility
//    NB_INTERNODE ,                        number of internodes (lateral axes)
    AUXILIARY                            // auxiliary variable (smoothing/piecewise linear function)
  };

  // Discrete parametric distributions

  const int MAX_INF_BOUND = 10000;       // maximum lower bound
  const int MAX_DIFF_BOUND = 10000;      // maximum difference between lower and upper bounds
  const double MAX_MEAN = 10000.;        // maximum mean

  const double B_PROBABILITY = 0.8;      // threshold for using the backward computation of the binomial probability mass function
  const double B_THRESHOLD = 1000.;      // threshold for using the computation in log of the binomial probability mass function
  const double P_THRESHOLD = 90.;        // threshold for using the computation in log of the Poisson probability mass function
  const double NB_THRESHOLD = 500.;      // threshold for using the computation in log of the negative binomial probability mass function

  const double SAMPLE_NB_VALUE_COEFF = 5.;  // factor for deducing the number of possible values of a distribution
                                            // from the number of possible values of a frequency distribution
  const int INF_BOUND_MARGIN = 5;        // range of values for the lower bound
  const int SUP_BOUND_MARGIN = 3;        // range of values for the upper bound
  const double POISSON_RATIO = 0.7;      // minimum mean/variance ratio for estimating a Poisson distribution
  const double POISSON_RANGE = 0.1;      // range for selecting a Poisson distribution by time scaling of
                                         // another discrete distribution
  const double NB_VALUE_COEFF = 2.;      // factor for deducing the number of values of a distribution from
                                         // the number of values of an initial distribution

  const int MIN_RANGE = 10;              // minimum interval of values for applying the rejection sampling method
  const double MAX_SURFACE = 3.;         // maximum surface for applying the rejection sampling method
  const int DIST_NB_ELEMENT = 1000000;   // maximum sample size for simulation

  // Renewal process estimated on the basis of time interval data

  const int NB_COMPLETE_INTERVAL = 3;    // minimum number of complete time intervals
  const double RENEWAL_LIKELIHOOD_DIFF = 1.e-5;  // threshold for stopping the EM iterations
  const int RENEWAL_NB_ITER = 10000;     // maximum number of EM iterations
  const double RENEWAL_DIFFERENCE_WEIGHT = 0.5;  // default penalty weight (1st- or 2nd-order difference cases)
  const double RENEWAL_ENTROPY_WEIGHT = 0.05;  // default penalty weight (entropy case)
  const double MAX_VALUE_COEFF = 10.;    // coefficient for deducing the maximum value of an inter-event distribution

  // Continuous parametric distributions

  const double CONTINUOUS_POSITIVE_INF_BOUND = 1.e-12; // inf bound of positive continuous distribution support (bug boost C++)

  const double GAMMA_TAIL = 1.e-3;       // gamma distribution tail
  const int GAMMA_NB_STEP = 1000;        // number of steps for computing a discretized gamma distribution
  const int GAMMA_NB_SUB_STEP = 10;      // number of sub-steps for computing a discretized gamma distribution
//  const int GAMMA_MIN_MEAN = 0.1;        // minimum mean of the gamma distribution
  const double GAMMA_INVERSE_SAMPLE_SIZE_FACTOR = 5.;  // factor for the gamma distribution corrected moment estimator
  const double GAMMA_MIN_SHAPE_PARAMETER = 0.1;  // minimum shape parameter (gamma distribution)
  const double GAMMA_DEFAULT_SCALE_PARAMETER = 1;  // default scale parameter (gamma distribution)
  const double GAMMA_ZERO_FREQUENCY_THRESHOLD = 0.999;  // threshold on the zero relative frequency (gamma distribution estimation)
  const double GAMMA_SHAPE_PARAMETER_THRESHOLD = 3.;  // threshold on the shape parameter (gamma distribution estimation)
  const double GAMMA_FREQUENCY_THRESHOLD = 100.;  // threshold on the frequency (gamma distribution estimation)
  const double GAMMA_ITERATION_FACTOR = 0.5;  // factor (gamma distribution estimation)
  const int GAMMA_MAX_NB_ITERATION = 5;  // maximum number of iterations (gamma distribution estimation)
//  const double GAMMA_VARIATION_COEFF_THRESHOLD = 1.e-2;  threshold on le coefficient de variation (gamma distribution estimation)

  const double INVERSE_GAUSSIAN_TAIL = 1.e-3;   // inverse Gaussian distribution tail
  const int INVERSE_GAUSSIAN_NB_STEP = 1000;    // number of steps for computing a discretized inverse Gaussian distribution
  const int INVERSE_GAUSSIAN_NB_SUB_STEP = 10;  // number of steps for computing a discretized inverse Gaussian distribution

  const double GAUSSIAN_TAIL = 5.e-4;    // Gaussian distribution tail
  const int GAUSSIAN_NB_STEP = 1000;     // number of steps for computing a discretized Gaussian distribution
  const int GAUSSIAN_NB_SUB_STEP = 10;   // number of steps for computing a discretized Gaussian distribution
  const double GAUSSIAN_MIN_VARIATION_COEFF = 1.e-3;  // minimum coefficient of variation (Gaussian distribution estimation)

  const int VON_MISES_NB_STEP = 3600;    // number of steps for computing a discretized von Mises distribution
  const int VON_MISES_NB_SUB_STEP = 10;  // number of steps for computing a discretized von Mises distribution
//  const double VON_MISES_CONCENTRATION_THRESHOLD = 10.;  threshold on the concentration parameter for applying
//                                                         the Gaussian approximation for computing a von Mises distribution

  const int CHI2_FREQUENCY = 2;          // minimum theoretical sample size for a Chi2 goodness of fit test

  const int MARGINAL_DISTRIBUTION_MAX_VALUE = 25000;  // maximum value for the building of a marginal frequency distribution
  const int HISTOGRAM_FREQUENCY = 10;    // average frequency for defining the bin width of an histogram
  const double SKEWNESS_ROUNDNESS = 1.e-2;  // rounding on the coefficient of skewness

  const int NB_ERROR = 10;               // maximum number of recorded errors

  const int LINE_NB_CHARACTER = 100;     // number of characters per line for sequences

  const int ASCII_NB_VALUE = 15;         // maximum number of values (ASCII display)
  const int ASCII_SPACE = 2;             // number of empty spaces between 2 columns (ASCII display)
  const double ASCII_ROUNDNESS = 1.e-5;  // rounding on the cumulative distribution function
                                         // for bounding a distribution (ASCII display)
  const double SPREADSHEET_ROUNDNESS = 1.e-7;  // rounding on the cumulative distribution function
                                               // for bounding a distribution (spreadsheet output)
  const int DISPLAY_NB_INDIVIDUAL = 50;  // maximum number of displayed individuals

  const int PLOT_NB_DISTRIBUTION = 15;   // maximum number of plotted distributions
  const int PLOT_NB_HISTOGRAM = 15;      // maximum number of plotted histograms
  const double PLOT_ROUNDNESS = 1.e-5;   // rounding on the cumulative distribution function
                                         // for bounding a plotted distribution
  const double PLOT_SHIFT = 0.2;         // shift between 2 plotted frequency distributions
  const double PLOT_MAX_SHIFT = 0.5;     // maximum shift between the first and last plotted frequency distributions
  const int TIC_THRESHOLD = 10;          // minimum number of plotted graduations
  const double PLOT_MASS_THRESHOLD = 1.e-3;  // minimum value for plotting a zero mass after the last possible value
  const double YSCALE = 1.1;             // scale factor for y axis in plots
  const double PLOT_RANGE_RATIO = 4.;    // threshold for plotting from 0



/****************************************************************
 *
 *  Class definition
 */


  /// \brief Test of hypothesis

  class Test {

    friend std::ostream& operator<<(std::ostream &os , const Test &test)
    { return test.ascii_print(os); }

  public :

    test_distribution ident;  ///< identifier (Normal/Chi2/F/Student's t)
    bool one_side;          ///< one-sided/two-sided
    int df1;                ///< degrees of freedom (Chi2/F/Student's t)
    int df2;                ///< degrees of freedom (F)
    double value;           ///< statistic value
    double critical_probability;  ///< critical probability

    void copy(const Test &test);

    Test(test_distribution iident , bool ione_side = true);
    Test(test_distribution iident , bool ione_side , int idf1 , int idf2 , double ivalue);
    Test(test_distribution iident , bool ione_side , int idf1 , int idf2 , double ivalue ,
         double icritical_probability);
    Test(const Test &test)
    { copy(test); }
    Test& operator=(const Test &test);

    std::ostream& ascii_print(std::ostream &os , bool comment_flag = false ,
                              bool reference_flag = true) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool reference_flag = true) const;

    void standard_normal_critical_probability_computation();
    void standard_normal_value_computation();
    void chi2_critical_probability_computation();
    void chi2_value_computation();
    void F_critical_probability_computation();
    void F_value_computation();
    void t_critical_probability_computation();
    void t_value_computation();
  };


  /// \brief Class for error management

  class StatError {

    friend std::ostream& operator<<(std::ostream &os , const StatError &error)
    { return error.ascii_write(os , ERROR); }

  private :

    int nb_error;           ///< number of errors
    int max_nb_error;       ///< maximum  number of errors
    int *line;              ///< line indices
    int *column;            ///< column indices
    char **label;           ///< error messages

  public :

    StatError(int imax_nb_error = NB_ERROR);
    ~StatError();

    std::ostream& ascii_write(std::ostream &os , error_type type = ERROR) const;

    void init() { nb_error = 0; }
    void update(const char *ilabel , int iline = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , const char *correction ,
                           int iline = 0 , int icolumn = 0);
    void correction_update(const char *ilabel , int correction ,
                           int iline = 0 , int icolumn = 0);

    // class member access

    int get_nb_error() const { return nb_error; }
    int get_max_nb_error() const { return max_nb_error; }
  };


  /// \brief Abstract class defining a common interface

  class StatInterface {

  public :

    virtual std::ostream& line_write(std::ostream &os) const = 0;

    virtual std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const = 0;
    std::string ascii_write(bool exhaustive = false) const;
    virtual bool ascii_write(StatError &error , const std::string path ,
                             bool exhaustive = false) const = 0;
    virtual bool spreadsheet_write(StatError &error , const std::string path) const = 0;

    virtual bool plot_write(StatError &error , const char *prefix ,
                            const char *title = NULL) const = 0;

    virtual MultiPlotSet* get_plotable() const { return NULL; };

//    bool binary_write(StatError &error , const std::string path) const;
  };


  class FrequencyDistribution;
  class DiscreteParametricModel;

  /// \brief Discrete distribution

  class Distribution {

    friend std::ostream& operator<<(std::ostream& , const Distribution&);

  public :

    int nb_value;           ///< number of values from 0
    int alloc_nb_value;     ///< number of allocated values
    int offset;             ///< number of values of null probability from 0
    double max;             ///< maximum probability
    double complement;      ///< complementary probability (> 0. for improper distributions)
    double mean;            ///< mean
    double variance;        ///< variance
    int nb_parameter;       ///< number of free parameters
    double *mass;           ///< probability mass function
    double *cumul;          ///< cumulative distribution function

    void mass_copy(const Distribution &dist , int inb_value = I_DEFAULT);
    void equal_size_copy(const Distribution &dist);
    void init(int inb_value);
    void copy(const Distribution &dist , int ialloc_nb_value = I_DEFAULT);
    void normalization_copy(const Distribution &dist);

    Distribution(int inb_value = 0);
    Distribution(int inb_value , double *imass);
    Distribution(const Distribution &dist , double scaling_coeff);
    Distribution(const FrequencyDistribution &histo);
    Distribution(const Distribution &dist ,
                 distribution_transformation transform = DISTRIBUTION_COPY ,
                 int ialloc_nb_value = I_DEFAULT);
    virtual ~Distribution();
    Distribution& operator=(const Distribution&);
    bool operator==(const Distribution&) const;
    bool operator!=(const Distribution &dist) const { return !(*this == dist); }

    std::ostream& ascii_characteristic_print(std::ostream &os , bool shape = false ,
                                             bool comment_flag = false) const;
    std::ostream& ascii_print(std::ostream &os , bool comment_flag ,
                              bool cumul_flag , bool nb_value_flag ,
                              const FrequencyDistribution *histo = NULL) const;

    std::ostream& ascii_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                              double *dist_scale , bool comment_flag , bool cumul_flag ,
                              const FrequencyDistribution *histo = NULL ,
                              bool mass_first = false) const;

    std::ostream& print(std::ostream &os) const;

    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;

    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag ,
                                    bool concentration_flag , bool nb_value_flag ,
                                    const FrequencyDistribution *histo = NULL) const;

    std::ostream& spreadsheet_print(std::ostream &os , int nb_dist , const Distribution **dist ,
                                    double *dist_scale , bool cumul_flag ,
                                    const FrequencyDistribution *histo = NULL ,
                                    bool mass_first = false) const;

    int plot_nb_value_computation(const FrequencyDistribution *histo = NULL) const;
    bool plot_print(const char *path , double *concentration , double scale) const;
    bool plot_print(const char *path , const FrequencyDistribution *histo = NULL) const;
    virtual std::ostream& plot_title_print(std::ostream &os) const
    { return os; }
    bool survival_plot_print(const char *path , double *survivor) const;

    bool plot_write(StatError &error , const char *prefix , int nb_dist ,
                    const Distribution **idist , const char *title) const;

    void plotable_mass_write(SinglePlot &plot , double scale = 1.) const;
    void plotable_cumul_write(SinglePlot &plot) const;
    void plotable_cumul_matching_write(SinglePlot &plot , const Distribution &reference_dist) const;
    void plotable_concentration_write(SinglePlot &plot) const;
    void plotable_survivor_write(SinglePlot &plot) const;

    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_distributions(StatError &error , int nb_dist ,
                                             const Distribution **idist) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(StatError &error , const std::string path) const;
    bool survival_spreadsheet_write(StatError &error , const std::string path) const;
    bool survival_plot_write(StatError &error , const char *prefix ,
                             const char *title = NULL) const;
    MultiPlotSet* survival_get_plotable(StatError &error) const;

    void max_computation();
    double mode_computation() const;
    void mean_computation();
    double quantile_computation(double icumul = 0.5) const;
    void variance_computation();
    void nb_value_computation();
    void offset_computation();
    double concentration_computation() const;

    double mean_absolute_deviation_computation(double location) const;
    double skewness_computation() const;
    double kurtosis_computation() const;
    double information_computation() const;

    double first_difference_norm_computation() const;
    double second_difference_norm_computation() const;

    void cumul_computation();
    double* survivor_function_computation() const;
    double* concentration_function_computation() const;

    double overlap_distance_computation(const Distribution &dist) const;

    void log_computation();

    double survivor_likelihood_computation(const FrequencyDistribution &histo) const;
    double chi2_value_computation(const FrequencyDistribution &histo) const;
    void chi2_degree_of_freedom(const FrequencyDistribution &histo , Test &test) const;

    void penalty_computation(double weight , penalty_type pen_type , double *penalty ,
                             side_effect outside) const;

    double likelihood_computation(const Reestimation<int> &histo) const
    { return histo.likelihood_computation(*this); }
    double likelihood_computation(const Reestimation<double> &histo) const
    { return histo.likelihood_computation(*this); }
    void chi2_fit(const FrequencyDistribution &histo , Test &test) const;

    void convolution(Distribution &dist1 , Distribution &dist2 ,
                     int inb_value = I_DEFAULT);
    int simulation() const;

    DiscreteParametricModel* truncate(StatError &error , int imax_value) const;
  };


  bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                  double *scale , int *dist_nb_value , int nb_histo ,
                  const FrequencyDistribution **histo);


  class Forward;

  /// \brief Discrete parametric distribution

  class DiscreteParametric : public Distribution {

    friend std::ostream& operator<<(std::ostream& , const DiscreteParametric&);

  public :

    discrete_parametric ident;  ///< identifier
    union {
      int inf_bound;        ///< lower bound
      int no_segment;       ///< number of segments (prior segment length distribution)
    };
    union {
      int sup_bound;        ///< upper bound (binomial, uniform)
      int sequence_length;  ///< sequence length (prior segment length distribution)
    };
    double parameter;       ///< parameter (Poisson, negative binomial, Poisson geometric)
    double probability;     ///< probability of success (binomial, negative binomial, Poisson geometric)

    void init(int iinf_bound , int isup_bound , double iparameter , double iprobability);
    void init(discrete_parametric iident , int iinf_bound , int isup_bound ,
              double iparameter , double iprobability);
    void copy(const DiscreteParametric &dist);

    DiscreteParametric(int inb_value = 0 , discrete_parametric iident = CATEGORICAL ,
                       int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
                       double iparameter = D_DEFAULT , double iprobability = D_DEFAULT);
    DiscreteParametric(discrete_parametric iident , int iinf_bound , int isup_bound ,
                       double iparameter , double iprobability ,
                       double cumul_threshold = CUMUL_THRESHOLD);
    DiscreteParametric(const Distribution &dist , int ialloc_nb_value = I_DEFAULT);
    DiscreteParametric(const Distribution &dist , double scaling_coeff);
    DiscreteParametric(const DiscreteParametric &dist , double scaling_coeff);
    DiscreteParametric(const FrequencyDistribution &histo);
    DiscreteParametric(const DiscreteParametric &dist ,
                       distribution_transformation transform = DISTRIBUTION_COPY ,
                       int ialloc_nb_value = I_DEFAULT);
    DiscreteParametric& operator=(const DiscreteParametric &dist);

    static DiscreteParametric* parsing(StatError &error , std::ifstream &in_file , int &line ,
                                       discrete_parametric last_ident = NEGATIVE_BINOMIAL ,
                                       double cumul_threshold = CUMUL_THRESHOLD ,
                                       int min_inf_bound = 0);

    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& ascii_parametric_characteristic_print(std::ostream &os , bool shape = false ,
                                                        bool comment_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    std::ostream& spreadsheet_parametric_characteristic_print(std::ostream &os , bool shape = false) const;
    std::ostream& plot_title_print(std::ostream &os) const;

    static int nb_value_computation(discrete_parametric ident , int inf_bound , int sup_bound ,
                                    double parameter , double probability ,
                                    double cumul_threshold = CUMUL_THRESHOLD);

    int nb_parameter_computation();
    void nb_parameter_update();

    double parametric_mean_computation() const;
    double parametric_variance_computation() const;
    double parametric_skewness_computation() const;
    double parametric_kurtosis_computation() const;

    double sup_norm_distance_computation(const DiscreteParametric &dist) const;

    void binomial_computation(int inb_value , distribution_computation mode);
    void poisson_computation(int inb_value , double cumul_threshold ,
                             distribution_computation mode);
    void negative_binomial_computation(int inb_value , double cumul_threshold ,
                                       distribution_computation mode);
    void poisson_geometric_computation(int inb_value , double cumul_threshold);
    void uniform_computation();
    void prior_segment_length_computation();

    void computation(int min_nb_value = 1 ,
                     double cumul_threshold = CUMUL_THRESHOLD);
    int simulation() const;

    double renewal_likelihood_computation(const Forward &forward_dist ,
                                          const FrequencyDistribution &within ,
                                          const FrequencyDistribution &backward ,
                                          const FrequencyDistribution &forward ,
                                          const FrequencyDistribution *no_event) const;
    void expectation_step(const FrequencyDistribution &within ,  // sequence_analysis
                          const FrequencyDistribution &backward ,
                          const FrequencyDistribution &forward ,
                          const FrequencyDistribution *no_event ,
                          Reestimation<double> *inter_event_reestim ,
                          Reestimation<double> *length_bias_reestim , int iter) const;

    double state_occupancy_likelihood_computation(const FrequencyDistribution &sojourn_time ,
                                                  const FrequencyDistribution &final_run) const;
    double state_occupancy_likelihood_computation(const Forward &forward ,  // sequence_analysis
                                                  const FrequencyDistribution &sojourn_time ,
                                                  const FrequencyDistribution &initial_run ,
                                                  const FrequencyDistribution &final_run ,
                                                  const FrequencyDistribution &single_run) const;
    void expectation_step(const FrequencyDistribution &sojourn_time ,  // sequence_analysis
                          const FrequencyDistribution &final_run ,
                          Reestimation<double> *occupancy_reestim , int iter) const;
    void expectation_step(const FrequencyDistribution &sojourn_time ,
                          const FrequencyDistribution &initial_run ,
                          const FrequencyDistribution &final_run ,
                          const FrequencyDistribution &single_run ,
                          Reestimation<double> *occupancy_reestim ,
                          Reestimation<double> *length_bias_reestim ,
                          int iter , bool combination = false ,
                          duration_distribution_mean_estimator mean_estimator = COMPUTED) const;
  };


  /// \brief Forward recurrence or sojourn time distribution

  class Forward : public DiscreteParametric {

  public :

    Forward(int inb_value = 0 , discrete_parametric iident = CATEGORICAL ,
            int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT ,
            double iparameter = D_DEFAULT , double iprobability = D_DEFAULT)
    :DiscreteParametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Forward(const DiscreteParametric &dist , int ialloc_nb_value = I_DEFAULT)
    :DiscreteParametric(dist , DISTRIBUTION_COPY , ialloc_nb_value) { computation(dist); }
    Forward(const Forward &forward , int ialloc_nb_value = I_DEFAULT)
    :DiscreteParametric((DiscreteParametric&)forward , DISTRIBUTION_COPY , ialloc_nb_value) {}

    void computation(const DiscreteParametric &dist);
  };


  class DiscreteDistributionData;
  class ContinuousParametric;
  class Convolution;
  class Compound;
  class DiscreteMixture;

  /// \brief Frequency distribution

  class FrequencyDistribution : public Reestimation<int> {

  public :

    FrequencyDistribution(int inb_value = 0)
    :Reestimation<int>(inb_value) {}
    FrequencyDistribution(const Distribution &dist)
    :Reestimation<int>(dist.nb_value) {}
    FrequencyDistribution(int inb_element , int *pelement);
    FrequencyDistribution(int nb_histo , const FrequencyDistribution **histo)
    :Reestimation<int>(nb_histo , (const Reestimation<int>**)histo) {}
    FrequencyDistribution(const FrequencyDistribution &histo ,
                          frequency_distribution_transformation transform ,
                          int param , rounding mode = FLOOR);

    bool operator==(const FrequencyDistribution&) const;
    bool operator!=(const FrequencyDistribution &histo) const { return !(*this == histo); }

    std::ostream& ascii_print(std::ostream &os , int comment_flag = false , bool cumul_flag = false) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_characteristic_print(std::ostream &os , bool shape = false) const;
    std::ostream& spreadsheet_circular_characteristic_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag = false ,
                                    bool concentration_flag = false) const;
    bool plot_print(const char *path , int nb_histo = 0 ,
                    const FrequencyDistribution **histo = NULL) const;
    bool plot_print(const char *path , double *cumul , double *concentration ,
                    double shift = 0.) const;
    bool survival_plot_print(const char *path , double *survivor) const;
    std::ostream& plot_title_print(std::ostream &os) const
    { return os; }

    void plotable_frequency_write(SinglePlot &plot) const;
    void plotable_mass_write(SinglePlot &plot) const;
    void plotable_cumul_write(SinglePlot &plot , double *icumul = NULL ,
                              double scale = D_DEFAULT) const;
    void plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,
                                       int  reference_nb_value , double *reference_cumul ,
                                       double *icumul = NULL) const;
    void plotable_concentration_write(SinglePlot &plot , double *icumul = NULL ,
                                      double scale = D_DEFAULT) const;
    void plotable_survivor_write(SinglePlot &plot) const;

    double* cumul_computation(double scale = D_DEFAULT) const;
    double concentration_computation() const;

    double* survivor_function_computation(double scale = D_DEFAULT) const;
    double* concentration_function_computation(double scale = D_DEFAULT) const;

    Test* kruskal_wallis_test(int nb_histo , const FrequencyDistribution **ihisto) const;

    std::ostream& dissimilarity_ascii_write(std::ostream &os , int nb_histo ,
                                            const FrequencyDistribution **ihisto ,
                                            variable_type type , double **dissimilarity) const;
    bool dissimilarity_ascii_write(StatError &error , const std::string path ,
                                   int nb_histo , const FrequencyDistribution **ihisto ,
                                   variable_type type , double **dissimilarity) const;
    bool dissimilarity_spreadsheet_write(StatError &error , const std::string path ,
                                         int nb_histo , const FrequencyDistribution **ihisto ,
                                         variable_type type , double **dissimilarity) const;

    void update(const Reestimation<double> *reestim , int inb_element);
    FrequencyDistribution* frequency_scale(int inb_element) const;
    double* rank_computation() const;
    int cumulative_distribution_function_computation(double **cdf) const;
    int min_interval_computation() const;

    DiscreteParametric* parametric_estimation(discrete_parametric ident , int min_inf_bound = 0 , bool flag = true ,
                                              double cumul_threshold = CUMUL_THRESHOLD) const;

    double likelihood_computation(const ContinuousParametric &dist ,
                                  int min_interval = I_DEFAULT) const;

    DiscreteDistributionData* merge(int nb_sample , std::vector<FrequencyDistribution> ihisto) const;

    void shift(const FrequencyDistribution &histo , int shift_param);
    void cluster(const FrequencyDistribution &histo , int step , rounding mode);

    DiscreteDistributionData* shift(StatError &error , int shift_param) const;
    DiscreteDistributionData* cluster(StatError &error , int step , rounding mode = FLOOR) const;
    DiscreteDistributionData* cluster(StatError &error , double ratio , bool display) const;
    DiscreteDistributionData* cluster(StatError &error , int nb_class , int *ilimit) const;
    DiscreteDistributionData* cluster(StatError &error , int nb_class , std::vector<int> ilimit) const;
    DiscreteDistributionData* transcode(StatError &error , int *category) const;
    DiscreteDistributionData* transcode(StatError &error , std::vector<int> category) const;
    DiscreteDistributionData* value_select(StatError &error , int min_value ,
                                           int max_value , bool keep = true) const;

    bool ascii_write(StatError &error , const std::string path) const;

    bool plot_write(StatError &error , const char *prefix , int nb_histo ,
                    const FrequencyDistribution **ihisto , const char *title) const;

    MultiPlotSet* get_plotable() const;
    MultiPlotSet* get_plotable_frequency_distributions(StatError &error , int nb_histo ,
                                                       const FrequencyDistribution **ihisto) const;

    std::ostream& survival_ascii_write(std::ostream &os) const;
    bool survival_ascii_write(StatError &error , const std::string path) const;
    bool survival_spreadsheet_write(StatError &error , const std::string path) const;
    bool survival_plot_write(StatError &error , const char *prefix ,
                             const char *title = NULL) const;
    MultiPlotSet* survival_get_plotable(StatError &error) const;

    bool comparison(StatError &error , bool display , int nb_histo ,
                    const FrequencyDistribution **ihisto , variable_type type ,
                    const std::string path = NULL , output_format format = ASCII) const;
    bool comparison(StatError &error , bool display , int nb_histo ,
                    const std::vector<FrequencyDistribution> ihisto , variable_type type ,
                    const std::string path = NULL , output_format format = ASCII) const;

    void F_comparison(bool display , const FrequencyDistribution &histo) const;
    void t_comparison(bool display , const FrequencyDistribution &histo) const;
    bool wilcoxon_mann_whitney_comparison(StatError &error , bool display ,
                                          const FrequencyDistribution &ihisto) const;

    DiscreteParametricModel* fit(StatError &error , const DiscreteParametric &idist) const;

    DiscreteParametricModel* parametric_estimation(StatError &error , discrete_parametric ident ,
                                                   int min_inf_bound = 0 , bool flag = true ,
                                                   double cumul_threshold = CUMUL_THRESHOLD) const;
    DiscreteParametricModel* type_parametric_estimation(StatError &error ,
                                                        int min_inf_bound = 0 , bool flag = true ,
                                                        double cumul_threshold = CUMUL_THRESHOLD) const;

    DiscreteMixture* discrete_mixture_estimation(StatError &error , const DiscreteMixture &imixt , bool *estimate ,
                                                 int min_inf_bound = 0 , bool mixt_flag = true ,
                                                 bool component_flag = true , double weight_step = 0.1) const;
    DiscreteMixture* discrete_mixture_estimation(StatError &error , const DiscreteMixture &imixt ,
                                                 int min_inf_bound = 0 , bool mixt_flag = true ,
                                                 bool component_flag = true , double weight_step = 0.1) const;
    DiscreteMixture* discrete_mixture_estimation(StatError &error , int nb_component , discrete_parametric *ident ,
                                                 int min_inf_bound = 0 , bool mixt_flag = true ,
                                                 bool component_flag = true , double weight_step = 0.1) const;
    DiscreteMixture* discrete_mixture_estimation(StatError &error , int nb_component , std::vector<discrete_parametric> ident ,
                                                 int min_inf_bound = 0 , bool mixt_flag = true ,
                                                 bool component_flag = true , double weight_step = 0.1) const;
    DiscreteMixture* discrete_mixture_estimation(StatError &error , bool display , int min_nb_component ,
                                                 int max_nb_component , discrete_parametric *ident , int min_inf_bound = 0 ,
                                                 bool mixt_flag = true , bool component_flag = true ,
                                                 model_selection_criterion criterion = BICc ,
                                                 double weight_step = 0.1) const;
    DiscreteMixture* discrete_mixture_estimation(StatError &error , bool display , int min_nb_component ,
                                                 int max_nb_component , std::vector<discrete_parametric> ident , int min_inf_bound = 0 ,
                                                 bool mixt_flag = true , bool component_flag = true ,
                                                 model_selection_criterion criterion = BICc ,
                                                 double weight_step = 0.1) const;

    Convolution* convolution_estimation(StatError &error , bool display , const DiscreteParametric &known_dist ,
                                        const DiscreteParametric &unknown_dist , estimation_criterion estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        penalty_type pen_type = SECOND_DIFFERENCE , side_effect outside = ZERO) const;
    Convolution* convolution_estimation(StatError &error , bool display , const DiscreteParametric &known_dist ,
                                        int min_inf_bound , estimation_criterion estimator = LIKELIHOOD ,
                                        int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                        penalty_type pen_type = SECOND_DIFFERENCE , side_effect outside = ZERO) const;

    Compound* compound_estimation(StatError &error , bool display , const DiscreteParametric &sum_dist ,
                                  const DiscreteParametric &dist , compound_distribution type ,
                                  estimation_criterion estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                  double weight = D_DEFAULT , penalty_type pen_type = SECOND_DIFFERENCE ,
                                  side_effect outside = ZERO) const;
    Compound* compound_estimation(StatError &error , bool display , const DiscreteParametric &known_dist ,
                                  compound_distribution type , int min_inf_bound = 0 ,
                                  estimation_criterion estimator = LIKELIHOOD ,
                                  int nb_iter = I_DEFAULT , double weight = D_DEFAULT ,
                                  penalty_type pen_type = SECOND_DIFFERENCE , side_effect outside = ZERO) const;

    DiscreteParametricModel* estimation(StatError &error , bool display ,
                                        const FrequencyDistribution &backward ,
                                        const FrequencyDistribution &forward ,
                                        const FrequencyDistribution *no_event ,
                                        const DiscreteParametric &iinter_event ,
                                        estimation_criterion estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                        duration_distribution_mean_estimator mean_estimator = COMPUTED ,
                                        double weight = D_DEFAULT , penalty_type pen_type = SECOND_DIFFERENCE ,
                                        side_effect outside = ZERO , double iinter_event_mean = D_DEFAULT) const;
    DiscreteParametricModel* estimation(StatError &error , bool display ,
                                        const FrequencyDistribution &backward ,
                                        const FrequencyDistribution &forward ,
                                        const FrequencyDistribution *no_event ,
                                        estimation_criterion estimator = LIKELIHOOD , int nb_iter = I_DEFAULT ,
                                        duration_distribution_mean_estimator mean_estimator = COMPUTED ,
                                        double weight = D_DEFAULT , penalty_type pen_type = SECOND_DIFFERENCE ,
                                        side_effect outside = ZERO) const;
  };


  bool plot_print(const char *path , int nb_dist , const Distribution **dist ,
                  double *scale , int *dist_nb_value , int nb_histo ,
                  const FrequencyDistribution **histo);


  class Histogram;

  /// \brief Continuous parametric distribution

  class ContinuousParametric {

//    friend std::ostream& operator<<(std::ostream &os , const ContinuousParametric &dist)
//    { return dist.ascii_print(os); }

  public :

    continuous_parametric ident;  ///< identifier
    union {
      double shape;         ///< shape parameter (gamma, zero-inflated gamma)
      double location;      ///< mean (Gaussian, inverse Gaussian, autoregressive model), mean direction (von Mises),
      double intercept;     ///< for linear model
    };
    union {
      double scale;         ///< scale parameter (gamma, inverse Gaussian, zero-inflated gamma)
      double dispersion;    ///< standard deviation (Gaussian, linear model, autoregressive model), concentration (von Mises),
    };
    union {
      double zero_probability;  ///< zero probability (zero-inflated gamma)
      double slope;         ///< for linear models
      double autoregressive_coeff;  ///< for autoregressive models
    };
    double min_value;       ///< minimum value
    double max_value;       ///< maximum value
    angle_unit unit;        ///< unit (degree/radian - von Mises)
    double slope_standard_deviation;  ///< for linear models
    double sample_size;     ///< for linear and autoregressive models
    union {
      double correlation;     ///< for linear models
      double determination_coeff;  ///< for autoregressive models
    };
    double *cumul;          ///< cumulative distribution function (von Mises)

    void copy(const ContinuousParametric &dist);

    ContinuousParametric(continuous_parametric iident = GAUSSIAN , double ilocation = D_INF ,
                         double idispersion = D_DEFAULT , double izero_probability = D_DEFAULT ,
                         angle_unit iunit = DEGREE);
    ContinuousParametric(const ContinuousParametric &dist)
    { copy(dist); }
    ~ContinuousParametric();
    ContinuousParametric& operator=(const ContinuousParametric&);

    static ContinuousParametric* parsing(StatError &error , std::ifstream &in_file ,
                                         int &line , continuous_parametric last_ident);

    std::ostream& ascii_parameter_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& ascii_characteristic_print(std::ostream &os ,
                                             bool file_flag = false) const;
    std::ostream& ascii_print(std::ostream &os , bool file_flag = false ,
                              bool cumul_flag = false , const Histogram *histo1 = NULL ,
                              const FrequencyDistribution *histo2 = NULL);
    std::ostream& spreadsheet_parameter_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os , bool cumul_flag = false ,
                                    const Histogram *histo1 = NULL ,
                                    const FrequencyDistribution *histo2 = NULL);
    std::ostream& spreadsheet_characteristic_print(std::ostream &os) const;
    std::ostream& plot_title_print(std::ostream &os) const;
    bool plot_print(const char *path , const Histogram *histo1 = NULL ,
                    const FrequencyDistribution *histo2 = NULL);
    bool q_q_plot_print(const char *path , int nb_value ,
                        double **empirical_cdf) const;
    void plotable_write(SinglePlot &plot , const Histogram *histo1 = NULL ,
                        const FrequencyDistribution *histo2 = NULL);
    void q_q_plotable_write(SinglePlot &plot , int nb_value ,
                            double **empirical_cdf) const;

    double** q_q_plot_computation(int nb_value , double **cdf) const;

    int nb_parameter_computation() const;

    double von_mises_mass_computation(double inf , double sup) const;
    double mass_computation(double inf , double sup) const;
    void von_mises_cumul_computation();

    double sup_norm_distance_computation(ContinuousParametric &dist);

    double likelihood_computation(const FrequencyDistribution &histo , int min_interval) const
    { return histo.likelihood_computation(*this , min_interval); }

    double simulation();
  };


  /// \brief Histogram

  class Histogram {

  public :

    int nb_element;         ///< sample size
    int nb_bin;             ///< number of bins
    double bin_width;       ///< constant bin width
    int max;                ///< maximum frequency
    int *frequency;         ///< frequency for each bin
    int type;               ///< variable type (INT_VALUE/REAL_VALUE)
    double min_value;       ///< minimum value
    double max_value;       ///< maximum value

    void copy(const Histogram &histo);

    Histogram(int inb_bin = 0 , bool init_flag = true);
    Histogram(const FrequencyDistribution &histo);
    Histogram(const Histogram &histo)
    { copy(histo); }
    ~Histogram();
    Histogram& operator=(const Histogram &histo);

    std::ostream& ascii_print(std::ostream &os , bool comment_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    bool plot_print(const char *path) const;
    void plotable_write(SinglePlot &plot) const;

    void max_computation();
    double* cumul_computation() const;
  };


  int column_width(int);
  int column_width(int min_value , int max_value);
  int column_width(double min_value , double max_value);
  int column_width(int nb_value , const double *value , double scale = 1.);
  char* label(const char *file_name);

  void cumul_computation(int nb_value , const double *pmass , double *pcumul);
  int cumul_method(int nb_value , const double *cumul , double scale = 1.);


};  // namespace stat_tool


#include "reestimation.hpp"



#endif
