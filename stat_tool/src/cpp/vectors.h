/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: vectors.h 18025 2015-04-23 07:09:13Z guedon $
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



namespace stat_tool {



/****************************************************************
 *
 *  Constantes :
 */


  const int VECTOR_NB_VARIABLE = 60;     // nombre maximum de variables
  const int DISTANCE_NB_VECTOR = 2000;   // nombre maximum de vecteurs pour le calcul
                                         // d'une matrice des distances
  const int CONTINGENCY_NB_VALUE = 100;  // nombre maximum de valeurs pour le calcul
                                         // d'un tableau de contingence
  const int DISPLAY_CONTINGENCY_NB_VALUE = 20;  // nombre maximum de valeurs pour l'affichage
                                                // d'un tableau de contingence
  const int VARIANCE_ANALYSIS_NB_VALUE = 100;  // nombre maximum de niveaux pour l'analyse de variance
  const int DISPLAY_CONDITIONAL_NB_VALUE = 100;  // nombre maximum de valeurs pour l'affichage
                                                 // des lois conditionnelles
  const int PLOT_NB_VALUE = 30;          // seuil pour l'ecriture des frequences (sortie Gnuplot)

  const int NB_SYMBOL = 50;              // nombre maximum de symboles

  const int MIN_NB_ASSIGNMENT = 1;   // nombre d'affectations des individus 1ere iteration de l'algorithme MCEM
  const int MAX_NB_ASSIGNMENT = 10;  // nombre d'affectations des individus  maximum pour l'algorithme MCEM
  const double NB_ASSIGNMENT_PARAMETER = 1.;  // parametre nombre d'affectations des individus pour l'algorithme MCEM

  enum {
    INDEPENDENT ,
    CONVOLUTION_FACTOR ,
    SCALING_FACTOR
  };

  enum {
    ABSOLUTE_VALUE ,
    QUADRATIC
  };



/****************************************************************
 *
 *  Definition des classes :
 */


  class ContinuousParametricProcess;
  class DistanceMatrix;
  class MultivariateMixture;
  class Mixture;
  class Regression;
  class VectorDistance;


  class Vectors : public StatInterface {  // vecteurs

    friend class Regression;
    friend class MultivariateMixture;
    friend class Mixture;

    friend Vectors* vectors_ascii_read(StatError &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const Vectors &vec)
    { return vec.ascii_write(os); }

  protected :

    int nb_vector;          // nombre de vecteurs
    int *identifier;        // identificateurs des vecteurs
    int nb_variable;        // nombre de variables
    int *type;              // type de chaque variable (INT_VALUE/REAL_VALUE)
    double *min_value;      // valeurs minimums
    double *max_value;      // valeurs maximums
    double *min_interval;   // intervalles minimums entre 2 valeurs
    FrequencyDistribution **marginal_distribution;  // lois marginales empiriques
    Histogram **marginal_histogram;  // histogrammes marginaux
    double *mean;           // vecteur moyenne
    double **covariance;    // matrice de variance-covariance
    int **int_vector;       // vecteurs, variables entieres [individual index][variable]
    double **real_vector;   // vecteurs, variables reelles  [individual index][variable]

    void init(int inb_vector , int *iidentifier , int inb_variable , int *itype , bool init_flag);
    void copy(const Vectors &vec);
    void add_state_variable(const Vectors &vec);
    void remove();

    void build_real_vector(int variable = I_DEFAULT);

    void transcode(const Vectors &vec , int variable , int min_symbol ,
                   int max_symbol , int *symbol);
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
    void build_marginal_histogram(int variable , double step = D_DEFAULT ,
                                  double imin_value = D_INF);
    int* order_computation(int variable) const;

    void mean_computation(int variable);
    void variance_computation(int variable);
    void covariance_computation(int variable = I_DEFAULT);

    double** correlation_computation() const;

    double** spearman_rank_correlation_computation() const;
    double** kendall_rank_correlation_computation() const;

    std::ostream& rank_correlation_ascii_write(std::ostream &os , int correlation_type ,
                                               double **correlation) const;
    bool rank_correlation_ascii_write(StatError &error , const char *path , int correlation_type ,
                                      double **correlation) const;

    int** joint_frequency_computation(int variable1 , int variable2) const;

    std::ostream& contingency_table_ascii_write(std::ostream &os , int variable1 , int variable2 ,
                                                int **frequency , double **deviation ,
                                                double **chi2_contribution , const Test &test ,
                                                bool file_flag = false) const;
    bool contingency_table_ascii_write(StatError &error , const char *path , int variable1 ,
                                       int variable2 , int **frequency , double **deviation ,
                                       double **chi2_contribution , const Test &test) const;
    bool contingency_table_spreadsheet_write(StatError &error , const char *path , int variable1 ,
                                             int variable2 , int **frequency , double **deviation ,
                                             double **chi2_contribution , const Test &test) const;

    std::ostream& variance_analysis_ascii_write(std::ostream &os , int type , const Vectors **value_vec ,
                                                bool exhaustive = false) const;
    bool variance_analysis_ascii_write(StatError &error , const char *path , int response_type ,
                                       const Vectors **value_vec , bool exhaustive = false) const;
    bool variance_analysis_spreadsheet_write(StatError &error , const char *path ,
                                             int response_type , const Vectors **value_vec) const;

    template <typename Type>
    void gamma_estimation(Type **component_vector_count , int variable ,
                          ContinuousParametricProcess *process , int iter) const;
    template <typename Type>
    void tied_gamma_estimation(Type **component_vector_count , int variable ,
                               ContinuousParametricProcess *process ,
                               int variance_factor , int iter) const;
    template <typename Type>
    void gaussian_estimation(Type **component_vector_count , int variable ,
                             ContinuousParametricProcess *process) const;
    template <typename Type>
    void tied_gaussian_estimation(Type **component_vector_count , int variable ,
                                  ContinuousParametricProcess *process ,
                                  int variance_factor) const;
    template <typename Type>
    void von_mises_estimation(Type **component_vector_count , int variable ,
                              ContinuousParametricProcess *process) const;

  public :

    Vectors();
    Vectors (int inb_vector , int *iidentifier , int inb_variable , int *itype ,
             bool init_flag = false)
    { init(inb_vector , iidentifier , inb_variable , itype , init_flag); }
    Vectors(int inb_vector , int *iidentifier , int inb_variable , int **iint_vector);  // interface AML
    Vectors(int inb_vector , int *iidentifier , int inb_variable , double **ireal_vector);  // interface AML
    Vectors(int inb_vector , int *iidentifier , int inb_variable , int *itype ,
            int **iint_vector , double **ireal_vector , bool variable_index = true);
    Vectors(const Vectors &vec , int variable , int itype);
    Vectors(const Vectors &vec , int inb_vector , int *index);
//    Vectors(const Vectors &vec , char transform = 'c' , int itype = I_DEFAULT);
    Vectors(const Vectors &vec , char transform = 'c');
    ~Vectors();
    Vectors& operator=(const Vectors &vec);

    void min_interval_computation(int variable);

    DiscreteDistributionData* extract(StatError &error , int variable) const;

    bool check(StatError &error);

    Vectors* merge(StatError &error, int nb_sample , const Vectors **ivec) const;
    Vectors* shift(StatError &error , int variable , int shift_param) const;
    Vectors* shift(StatError &error , int variable , double shift_param) const;
    Vectors* thresholding(StatError &error , int variable , int threshold , int mode) const;
    Vectors* thresholding(StatError &error , int variable , double threshold , int mode) const;
    Vectors* transcode(StatError &error , int variable , int *symbol) const;
    Vectors* cluster(StatError &error , int variable , int step ,
                     int mode = FLOOR) const;
    Vectors* cluster(StatError &error , int variable , int inb_value ,
                     int *ilimit) const;
    Vectors* cluster(StatError &error , int variable , int nb_class ,
                     double *ilimit) const;
    Vectors* scaling(StatError &error , int variable , int scaling_coeff) const;
    Vectors* scaling(StatError &error , int variable , double scaling_coeff) const;
    Vectors* round(StatError &error , int variable = I_DEFAULT ,
                   int mode = ROUND) const;

    Vectors* value_select(StatError &error , std::ostream &os , int variable ,
                          int imin_value , int imax_value , bool keep = true) const;
    Vectors* value_select(StatError &error , std::ostream &os , int variable ,
                          double imin_value , double imax_value , bool keep = true) const;

    Vectors* select_individual(StatError &error , int inb_vector , int *iidentifier ,
                               bool keep = true) const;
    Vectors* select_variable(StatError &error , int inb_variable , int *ivariable ,
                             bool keep = true) const;
    Vectors* merge_variable(StatError &error , int nb_sample , const Vectors **ivec ,
                            int ref_sample = I_DEFAULT) const;

    std::ostream& line_write(std::ostream &os) const;

    virtual std::ostream& ascii_data_write(std::ostream &os , bool exhaustive = false) const;
    virtual bool ascii_data_write(StatError &error , const char *path , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    bool select_step(StatError &error , int variable , double step ,
                     double imin_value = D_INF);
    int cumulative_distribution_function_computation(int variable , double **cdf) const;

    double mean_absolute_deviation_computation(int variable) const;
    double mean_absolute_difference_computation(int variable) const;
    double skewness_computation(int variable) const;
    double kurtosis_computation(int variable) const;

    double* mean_direction_computation(int variable , int unit) const;

    double spearman_rank_single_correlation_computation() const;
    double kendall_rank_single_correlation_computation() const;

    bool rank_correlation_computation(StatError &error , std::ostream &os ,
                                      int correlation_type , const char *path = NULL) const;

    DistanceMatrix* comparison(StatError &error , const VectorDistance &ivector_dist ,
                               bool standardization = true) const;

    bool contingency_table(StatError &error , std::ostream &os , int variable1 ,
                           int variable2 , const char *path = NULL , char format = 'a') const;

    bool variance_analysis(StatError &error , std::ostream &os , int class_variable ,
                           int response_variable , int response_type ,
                           const char *path = NULL , char format = 'a') const;

    Regression* linear_regression(StatError &error , int explanatory_variable ,
                                  int response_variable) const;
    Regression* moving_average(StatError &error , int explanatory_variable ,
                               int response_variable , int nb_point ,
                               double *filter , char algorithm = 'a') const;
    Regression* moving_average(StatError &error , int explanatory_variable ,
                               int response_variable , const Distribution &dist ,
                               char algorithm = 'a') const;
    Regression* nearest_neighbor_smoother(StatError &error , int explanatory_variable ,
                                          int response_variable , double span ,
                                          bool weighting = true) const;

    MultivariateMixture* mixture_estimation(StatError &error , std::ostream &os ,
                                            const MultivariateMixture &imixture,
                                            int nb_iter = I_DEFAULT ,
                                            bool *force_param = NULL) const;
    MultivariateMixture* mixture_estimation(StatError &error , std::ostream& os ,
                                            int nb_component , int nb_iter = I_DEFAULT ,
                                            bool *force_param = NULL) const;

    Mixture* mixture_estimation(StatError &error , std::ostream &os , const Mixture &imixt ,
                                bool known_component = false , bool common_dispersion = false ,
                                int variance_factor = INDEPENDENT , bool assignment = true ,
                                int nb_iter = I_DEFAULT) const;
    Mixture* mixture_estimation(StatError &error , ostream &os , int nb_component ,
                                int ident , double mean , double standard_deviation ,
                                bool tied_location = true , int variance_factor = SCALING_FACTOR ,
                                bool assignment = true , int nb_iter = I_DEFAULT) const;
    Mixture* mixture_stochastic_estimation(StatError &error , std::ostream &os , const Mixture &imixt ,
                                           bool known_component = false , bool common_dispersion = false ,
                                           int variance_factor = INDEPENDENT ,
                                           int min_nb_assignment = MIN_NB_ASSIGNMENT ,
                                           int max_nb_assignment = MAX_NB_ASSIGNMENT ,
                                           double parameter = NB_ASSIGNMENT_PARAMETER ,
                                           bool assignment = true , int nb_iter = I_DEFAULT) const;
    Mixture* mixture_stochastic_estimation(StatError &error , ostream &os , int nb_component ,
                                           int ident , double mean , double standard_deviation ,
                                           bool tied_location = true , int variance_factor = SCALING_FACTOR ,
                                           int min_nb_assignment = MIN_NB_ASSIGNMENT ,
                                           int max_nb_assignment = MAX_NB_ASSIGNMENT ,
                                           double parameter = NB_ASSIGNMENT_PARAMETER ,
                                           bool assignment = true , int nb_iter = I_DEFAULT) const;

    // acces membres de la classe

    int get_nb_vector() const { return nb_vector; }
    int get_identifier(int ivec) const { return identifier[ivec]; }
    int get_nb_variable() const { return nb_variable; }
    int get_type(int variable) const { return type[variable]; }
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


  Vectors* vectors_ascii_read(StatError &error , const char *path);



  class VectorDistance : public StatInterface {  // parametres de definition
                                                 // d'une distance entre vecteurs

    friend class Vectors;

    friend VectorDistance* vector_distance_ascii_read(StatError &error , const char *path);
    friend std::ostream& operator<<(std::ostream &os , const VectorDistance &param)
    { return param.ascii_write(os); }

  private :

    int nb_variable;        // nombre de variables
    int distance_type;      // type de distance (ABSOLUTE_VALUE/QUADRATIC)
    int *variable_type;     // type de chaque variable (SYMBOLIC/ORDINAL/NUMERIC/CIRCULAR)
    double *weight;         // poids de chaque variable
    double *dispersion;     // quantite pour la standardisation
    int *nb_value;          // nombre de valeurs par variable
    double ***symbol_distance;  // matrice des distances entre symboles
    int *period;            // periode (variable circulaire)

    void copy(const VectorDistance &param);
    void remove();

  public :

    VectorDistance();
    VectorDistance(int inb_variable , int *ivariable_type , double *iweight ,
                   int idistance_type = ABSOLUTE_VALUE);
    VectorDistance(int inb_variable , int idistance_type , int *ivariable_type ,
                   double *iweight , int *inb_value , double ***isymbol_distance ,
                   int *iperiod);
    VectorDistance(const VectorDistance &vector_dist)
    { copy(vector_dist); }
    ~VectorDistance();
    VectorDistance& operator=(const VectorDistance &vector_dist);

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;

    // fonctions pour la compatibilite avec la classe StatInterface

    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = NULL) const;

    double* max_symbol_distance_computation(int variable) const;

    void dispersion_update(int variable , double idispersion) const;
    void dispersion_computation(int variable , const FrequencyDistribution *marginal_distribution ,
                                double *rank = NULL) const;

    // acces membres de la classe

    int get_nb_variable() const { return nb_variable; }
    int get_distance_type() const { return distance_type; }
    int get_variable_type(int variable) const { return variable_type[variable]; }
    double get_weight(int variable) const { return weight[variable]; }
    double get_dispersion(int variable) const { return dispersion[variable]; }
    int get_nb_value(int variable) const { return nb_value[variable]; }
    double** get_symbol_distance(int variable) const { return symbol_distance[variable]; }
    double get_symbol_distance(int variable , int symbol1 , int symbol2) const
    { return symbol_distance[variable][symbol1][symbol2]; }
    int get_period(int variable) const { return period[variable]; }
  };


  VectorDistance* vector_distance_ascii_read(StatError &error , const char *path);

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
