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
 *       $Id$
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



#ifndef VECTORS_H
#define VECTORS_H


#include "stat_tools.h"
#include "distribution.h"
#include "distance_matrix.h"
#include "markovian.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const int VECTOR_NB_VARIABLE = 60;     // maximum number of variables
  const int DISTANCE_NB_VECTOR = 2000;   // maximum number of vectors for the computation
                                         // of a matrix of pairwise distances
  const int SUP_NORM_DISTANCE_NB_VECTOR = 10;   // minimum number of vectors for the computation
                                                // of the sup norm distance
  const int CONTINGENCY_NB_VALUE = 100;  // maximum number of categories for the computation
                                         // of a contingency table
  const int DISPLAY_CONTINGENCY_NB_VALUE = 20;  // maximum number of categories for the display of
                                                // a contingency table
  const int VARIANCE_ANALYSIS_NB_VALUE = 100;  // maximum number of levels for the analysis of variance
  const int DISPLAY_CONDITIONAL_NB_VALUE = 100;  // maximum number of values of the display of
                                                 // the conditional distributions
  const int PLOT_NB_VALUE = 30;          // threshold for the writing of frequencies (Gnuplot output)

  const int NB_CATEGORY = 50;            // maximum number of categories

  const int MIN_NB_ASSIGNMENT = 1;   // minimum number of assignments of individuals (1st iteration of the MCEM algorithm)
  const int MAX_NB_ASSIGNMENT = 10;  // maximum number of assignments of individuals (MCEM algorithm)
  const double NB_ASSIGNMENT_PARAMETER = 1.;  // parameter for the number of assignments of individuals (MCEM algorithm)

  enum vector_transformation {
    VECTOR_COPY ,
    ADD_COMPONENT_VARIABLE
  };

  enum threshold_direction {
    ABOVE ,
    BELOW
  };

  enum correlation_type {
    PEARSON ,
    SPEARMAN ,
    KENDALL ,
    SPEARMAN2
  };

  enum metric {
    ABSOLUTE_VALUE ,
    QUADRATIC
  };

  enum tying_rule {
    INDEPENDENT ,
    CONVOLUTION_FACTOR ,
    SCALING_FACTOR
  };

  enum moving_average_method {
    AVERAGING ,
    LEAST_SQUARES
  };



/****************************************************************
 *
 *  Class definition
 */


  class Mixture;
  class Regression;
  class VectorDistance;

  /// \brief Vectors with integer- and real-valued variables

  class Vectors : public StatInterface {

    friend class Regression;
    friend class Mixture;

    friend std::ostream& operator<<(std::ostream &os , const Vectors &vec)
    { return vec.ascii_write(os); }

  protected :

    int nb_vector;          ///< number of vectors
    int *identifier;        ///< vector identifiers
    int nb_variable;        ///< number of variables
    variable_nature *type;  ///< variable types (INT_VALUE/REAL_VALUE)
    double *min_value;      ///< minimum value for each variable
    double *max_value;      ///< maximum value for each variable
    double *min_interval;   ///< minimum interval between 2 values for each variable
    FrequencyDistribution **marginal_distribution;  ///< marginal frequency distributions
    Histogram **marginal_histogram;  ///< marginal histograms
    double *mean;           ///< mean vector
    double **covariance;    ///< variance-covariance matrix
    int **int_vector;       ///< vectors, integer-valued variables
    double **real_vector;   ///< vectors, real-valued variables

    void init(int inb_vector , int *iidentifier , int inb_variable ,
              variable_nature *itype , bool init_flag);
    void copy(const Vectors &vec);
    void add_state_variable(const Vectors &vec);
    void remove();

    void build_real_vector(int variable = I_DEFAULT);

    void transcode(const Vectors &vec , int variable , int min_category ,
                   int max_category , int *category);
    void cluster(const Vectors &vec , int variable , int nb_class , double *limit);
    void select_variable(const Vectors &vec , int *variable);
    Vectors* remove_variable_1() const;
    Vectors* select_variable(int explanatory_variable , int response_variable) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive , bool comment_flag) const;
    std::ostream& ascii_print(std::ostream &os , bool comment_flag ,
                              double *posterior_probability = NULL , double *entropy = NULL) const;
    bool plot_print(const char *path , double *standard_residual = NULL) const;
    void plotable_write(SinglePlot &plot , int variable1 , int variable2) const;
    void plotable_frequency_write(SinglePlot &plot , int variable1 , int variable2) const;

    void min_value_computation(int variable);
    void max_value_computation(int variable);

    void build_marginal_frequency_distribution(int variable);
    void build_marginal_histogram(int variable , double bin_width = D_DEFAULT ,
                                  double imin_value = D_INF);
    int* order_computation(int variable) const;

    void mean_computation(int variable);
    void variance_computation(int variable);
    void covariance_computation(int variable = I_DEFAULT);

    double** correlation_computation() const;

    double** spearman_rank_correlation_computation() const;
    double** kendall_rank_correlation_computation() const;

    std::ostream& rank_correlation_ascii_write(std::ostream &os , correlation_type correl_type ,
                                               double **correlation) const;
    bool rank_correlation_ascii_write(StatError &error , const std::string path , correlation_type correl_type ,
                                      double **correlation) const;

    int** joint_frequency_computation(int variable1 , int variable2) const;

    std::ostream& contingency_table_ascii_write(std::ostream &os , int variable1 , int variable2 ,
                                                int **frequency , double **deviation ,
                                                double **chi2_contribution , const Test &test ,
                                                bool file_flag = false) const;
    bool contingency_table_ascii_write(StatError &error , const std::string path , int variable1 ,
                                       int variable2 , int **frequency , double **deviation ,
                                       double **chi2_contribution , const Test &test) const;
    bool contingency_table_spreadsheet_write(StatError &error , const std::string path , int variable1 ,
                                             int variable2 , int **frequency , double **deviation ,
                                             double **chi2_contribution , const Test &test) const;

    std::ostream& variance_analysis_ascii_write(std::ostream &os , int type , const Vectors **value_vec ,
                                                bool exhaustive = false) const;
    bool variance_analysis_ascii_write(StatError &error , const std::string path , int response_type ,
                                       const Vectors **value_vec , bool exhaustive = false) const;
    bool variance_analysis_spreadsheet_write(StatError &error , const std::string path ,
                                             int response_type , const Vectors **value_vec) const;

    double likelihood_computation(int variable , const ContinuousParametric &dist) const;

    double gamma_estimation(int variable , ContinuousParametric *dist) const;
    double inverse_gaussian_estimation(int variable , ContinuousParametric *dist) const;
    double gaussian_estimation(int variable , ContinuousParametric *dist) const;
    double von_mises_estimation(int variable , ContinuousParametric *dist) const;

    double continuous_parametric_estimation(int variable , continuous_parametric ident) const;

    template <typename Type>
    void gamma_estimation(Type **component_vector_count , int variable ,
                          ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void tied_gamma_estimation(Type **component_vector_count , int variable ,
                               ContinuousParametricProcess *process ,
                               tying_rule variance_factor , int iter) const;
    template <typename Type>
    void inverse_gaussian_estimation(Type **component_vector_count , int variable ,
                                     ContinuousParametricProcess *process) const;
    template <typename Type>
    void tied_inverse_gaussian_estimation(Type **component_vector_count , int variable ,
                                          ContinuousParametricProcess *process ,
                                          tying_rule variance_factor) const;
    template <typename Type>
    void gaussian_estimation(Type **component_vector_count , int variable ,
                             ContinuousParametricProcess *process) const;
    template <typename Type>
    void tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                  ContinuousParametricProcess *process) const;
    template <typename Type>
    void tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                  ContinuousParametricProcess *process ,
                                  tying_rule variance_factor) const;
    template <typename Type>
    void von_mises_estimation(Type **component_vector_count , int variable ,
                              ContinuousParametricProcess *process) const;

  public :

    Vectors();
    Vectors (int inb_vector , int *iidentifier , int inb_variable , variable_nature *itype ,
             bool init_flag = false)
    { init(inb_vector , iidentifier , inb_variable , itype , init_flag); }
    Vectors(int inb_vector , int *iidentifier , int inb_variable , int **iint_vector);  // AML interface
    Vectors(int inb_vector , int *iidentifier , int inb_variable , double **ireal_vector);  // AML interface
    Vectors(int inb_vector , int *iidentifier , int inb_variable , variable_nature *itype ,
            int **iint_vector , double **ireal_vector , bool variable_index = true);
    Vectors(const Vectors &vec , int variable , variable_nature itype);
    Vectors(const Vectors &vec , int inb_vector , int *index);
//    Vectors(const Vectors &vec , vector_transformation transform = VECTOR_COPY , variable_nature itype = I_DEFAULT);
    Vectors(const Vectors &vec , vector_transformation transform = VECTOR_COPY);
    ~Vectors();
    Vectors& operator=(const Vectors &vec);

    void min_interval_computation(int variable);

    DiscreteDistributionData* extract(StatError &error , int variable) const;

    bool check(StatError &error);

    Vectors* merge(StatError &error, int nb_sample , const Vectors **ivec) const;
    Vectors* merge(StatError &error, int nb_sample , const std::vector<Vectors> ivec) const;

    Vectors* shift(StatError &error , int variable , int shift_param) const;
    Vectors* shift(StatError &error , int variable , double shift_param) const;
    Vectors* thresholding(StatError &error , int variable , int threshold ,
                          threshold_direction mode) const;
    Vectors* thresholding(StatError &error , int variable , double threshold ,
                          threshold_direction mode) const;
    Vectors* transcode(StatError &error , int variable , int *category) const;
    Vectors* transcode(StatError &error , int variable , std::vector<int> category) const;

    Vectors* cluster(StatError &error , int variable , int step ,
                     rounding mode = FLOOR) const;
    Vectors* cluster(StatError &error , int variable , int inb_value ,
                     int *ilimit) const;
    Vectors* cluster(StatError &error , int variable , int inb_value ,
                     std::vector<int> ilimit) const;
    Vectors* cluster(StatError &error , int variable , int nb_class ,
                     double *ilimit) const;
    Vectors* cluster(StatError &error , int variable , int nb_class ,
                     std::vector<double> ilimit) const;
    Vectors* scaling(StatError &error , int variable , int scaling_coeff) const;
    Vectors* scaling(StatError &error , int variable , double scaling_coeff) const;
    Vectors* round(StatError &error , int variable = I_DEFAULT ,
                   rounding mode = ROUND) const;
    Vectors* log_transform(StatError &error , int variable = I_DEFAULT ,
                           log_base base = NATURAL) const;

    Vectors* value_select(StatError &error , bool display , int variable ,
                          int imin_value , int imax_value , bool keep = true) const;
    Vectors* value_select(StatError &error , bool display , int variable ,
                          double imin_value , double imax_value , bool keep = true) const;

    Vectors* select_individual(StatError &error , int inb_vector , int *iidentifier ,
                               bool keep = true) const;
    Vectors* select_individual(StatError &error , int inb_vector , std::vector<int> iidentifier ,
                               bool keep = true) const;
    Vectors* select_variable(StatError &error , int inb_variable , int *ivariable ,
                             bool keep = true) const;
    Vectors* select_variable(StatError &error , int inb_variable , std::vector<int> ivariable ,
                             bool keep = true) const;
    Vectors* sum_variable(StatError &error , int nb_summed_variable , int *ivariable) const;
    Vectors* sum_variable(StatError &error , int nb_summed_variable , std::vector<int> ivariable) const;
    Vectors* merge_variable(StatError &error , int nb_sample , const Vectors **ivec ,
                            int ref_sample = I_DEFAULT) const;
    Vectors* merge_variable(StatError &error , int nb_sample , const std::vector<Vectors> ivec ,
                            int ref_sample = I_DEFAULT) const;

    static Vectors* ascii_read(StatError &error , const std::string path);

    std::ostream& line_write(std::ostream &os) const;

    virtual std::ostream& ascii_data_write(std::ostream &os , bool exhaustive = false) const;
    virtual bool ascii_data_write(StatError &error , const std::string path , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    bool select_bin_width(StatError &error , int variable , double bin_width ,
                          double imin_value = D_INF);
    int cumulative_distribution_function_computation(int variable , double **cdf) const;

    double mean_absolute_deviation_computation(int variable , double location) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable) const;
    double kurtosis_computation(int variable) const;

    double* mean_direction_computation(int variable , angle_unit unit) const;

    double spearman_rank_single_correlation_computation() const;
    double kendall_rank_single_correlation_computation() const;

    bool rank_correlation_computation(StatError &error , bool display ,
                                      correlation_type correl_type , const std::string path = "") const;

    DistanceMatrix* comparison(StatError &error , const VectorDistance &ivector_dist ,
                               bool standardization = true) const;

    bool contingency_table(StatError &error , bool display , int variable1 ,
                           int variable2 , const std::string path = "" , output_format format = ASCII) const;

    bool variance_analysis(StatError &error , bool display , int class_variable ,
                           int response_variable , int response_type ,
                           const std::string path = "" , output_format format = ASCII) const;

    bool sup_norm_distance(StatError &error , bool display , const Vectors &ivec) const;

    Regression* linear_regression(StatError &error , int explanatory_variable ,
                                  int response_variable) const;
    Regression* moving_average(StatError &error , int explanatory_variable ,
                               int response_variable , int nb_point , double *filter ,
                               moving_average_method algorithm = AVERAGING) const;
    Regression* moving_average(StatError &error , int explanatory_variable ,
                               int response_variable , const Distribution &dist ,
                               moving_average_method algorithm = AVERAGING) const;
    Regression* nearest_neighbor_smoother(StatError &error , int explanatory_variable ,
                                          int response_variable , double span ,
                                          bool weighting = true) const;

    Mixture* mixture_estimation(StatError &error , bool display , const Mixture &imixt ,
                                bool known_component = false , bool common_dispersion = false ,
                                tying_rule variance_factor = INDEPENDENT , bool assignment = true ,
                                int nb_iter = I_DEFAULT) const;
    Mixture* mixture_estimation(StatError &error , bool display , int nb_component ,
                                double offset , double mean , double standard_deviation ,
                                bool common_dispersion = true , bool assignment = true ,
                                int nb_iter = I_DEFAULT) const;
    Mixture* mixture_estimation(StatError &error , bool display , int nb_component ,
                                int ident , double mean , double standard_deviation ,
                                bool tied_location = true , tying_rule variance_factor = SCALING_FACTOR ,
                                bool assignment = true , int nb_iter = I_DEFAULT) const;
    Mixture* mixture_stochastic_estimation(StatError &error , bool display , const Mixture &imixt ,
                                           bool known_component = false , bool common_dispersion = false ,
                                           tying_rule variance_factor = INDEPENDENT ,
                                           int min_nb_assignment = MIN_NB_ASSIGNMENT ,
                                           int max_nb_assignment = MAX_NB_ASSIGNMENT ,
                                           double parameter = NB_ASSIGNMENT_PARAMETER ,
                                           bool assignment = true , int nb_iter = I_DEFAULT) const;
    Mixture* mixture_stochastic_estimation(StatError &error , bool display , int nb_component ,
                                           double offset , double mean , double standard_deviation ,
                                           bool common_dispersion = true ,
                                           int min_nb_assignment = MIN_NB_ASSIGNMENT ,
                                           int max_nb_assignment = MAX_NB_ASSIGNMENT ,
                                           double parameter = NB_ASSIGNMENT_PARAMETER ,
                                           bool assignment = true , int nb_iter = I_DEFAULT) const;
    Mixture* mixture_stochastic_estimation(StatError &error , bool display , int nb_component ,
                                           int ident , double mean , double standard_deviation ,
                                           bool tied_location = true , tying_rule variance_factor = SCALING_FACTOR ,
                                           int min_nb_assignment = MIN_NB_ASSIGNMENT ,
                                           int max_nb_assignment = MAX_NB_ASSIGNMENT ,
                                           double parameter = NB_ASSIGNMENT_PARAMETER ,
                                           bool assignment = true , int nb_iter = I_DEFAULT) const;

    // class member access

    int get_nb_vector() const { return nb_vector; }
    int get_identifier(int ivec) const { return identifier[ivec]; }
    int get_nb_variable() const { return nb_variable; }
    variable_nature get_type(int variable) const { return type[variable]; }
    double get_min_value(int variable) const { return min_value[variable]; }
    double get_max_value(int variable) const { return max_value[variable]; }
    FrequencyDistribution* get_marginal_distribution(int variable) const
    { return marginal_distribution[variable]; }
    Histogram* get_marginal_histogram(int variable) const
    { return marginal_histogram[variable]; }
    double get_mean(int variable) const { return mean[variable]; }
    double get_covariance(int variable1, int variable2) const
    { return covariance[variable1][variable2]; }
    int get_int_vector(int ivec , int variable) const
    { return int_vector[ivec][variable]; }
    double get_real_vector(int ivec , int variable) const
    { return real_vector[ivec][variable]; }
  };


  /// \brief Parameterization of a distance between vectors with heterogeneous variables

  class VectorDistance : public StatInterface {

    friend class Vectors;

    friend std::ostream& operator<<(std::ostream &os , const VectorDistance &param)
    { return param.ascii_write(os); }

  private :

    int nb_variable;        ///< number of variables
    metric distance_type;   ///< distance type (ABSOLUTE_VALUE/QUADRATIC)
    variable_type *var_type;  ///< variable types (NOMINAL/ORDINAL/NUMERIC/CIRCULAR)
    double *weight;         ///< weight of each variable
    double *dispersion;     ///< dispersion quantities for the standardization
    int *nb_value;          ///< number of values for each variable
    double ***category_distance;  ///< matrices of distances between categories (for categorical variables)
    int *period;            ///< periods (for circular variables)

    void copy(const VectorDistance &param);
    void remove();

  public :

    VectorDistance();
    VectorDistance(int inb_variable , variable_type *ivar_type , double *iweight ,
                   metric idistance_type = ABSOLUTE_VALUE);
    VectorDistance(int inb_variable , metric idistance_type , variable_type *ivar_type ,
                   double *iweight , int *inb_value , double ***icategory_distance ,
                   int *iperiod);
    VectorDistance(const VectorDistance &vector_dist)
    { copy(vector_dist); }
    ~VectorDistance();
    VectorDistance& operator=(const VectorDistance &vector_dist);

    static VectorDistance* ascii_read(StatError &error , const std::string path);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;

    // functions for the compatibility with the StatInterface class

    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;

    double* max_category_distance_computation(int variable) const;

    void dispersion_update(int variable , double idispersion) const;
    void dispersion_computation(int variable , const FrequencyDistribution *marginal_distribution ,
                                double *rank = NULL) const;

    // class member access

    int get_nb_variable() const { return nb_variable; }
    metric get_distance_type() const { return distance_type; }
    variable_type get_var_type(int variable) const { return var_type[variable]; }
    double get_weight(int variable) const { return weight[variable]; }
    double get_dispersion(int variable) const { return dispersion[variable]; }
    int get_nb_value(int variable) const { return nb_value[variable]; }
    double** get_category_distance(int variable) const { return category_distance[variable]; }
    double get_category_distance(int variable , int category1 , int category2) const
    { return category_distance[variable][category1][category2]; }
    int get_period(int variable) const { return period[variable]; }
  };

  bool identifier_checking(StatError &error , int nb_individual , int *identifier);
  bool selected_identifier_checking(StatError &error , int nb_individual , int *identifier ,
                                    int nb_selected_individual , int *selected_identifier ,
                                    const char *data_label);
  int* identifier_select(int nb_individual , int *identifier , int nb_selected_individual ,
                         int *selected_identifier , bool keep);
  int* select_variable(int nb_variable , int nb_selected_variable ,
                       int *selected_variable , bool keep);


};  // namespace stat_tool


#include "continuous_parametric_estimation.hpp"



#endif
