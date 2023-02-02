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



#ifndef MARKOVIAN_H
#define MARKOVIAN_H


#include "stat_tools.h"
#include "chain_reestimation.h"


namespace stat_tool {



/****************************************************************
 *
 *  Constants
 */


  const int NB_STATE = 100;              // maximum number of states of a Markov chain
  const int ORDER = 8;                   // maximum order of a Markov chain
  const double MIN_PROBABILITY = 1.e-5;  // minimum initial/transition/categorical observation probability
  const double THRESHOLDING_FACTOR = 0.8;  // factor for the thresholding of probabilities
  const int NB_PARAMETER = 100000;       // maximum number of parameters of a Markov chain
  const int NB_OUTPUT_PROCESS = 15;      // maximum number of observation processes
  const int NB_OUTPUT = 25;              // maximum number of observed categories per state (categorical case)
  const double OBSERVATION_THRESHOLD = 0.999;  // threshold on the cumulative distribution function for bounding
                                               // a discrete parametric observation distribution

  const double ACCESSIBILITY_THRESHOLD = 1.e-6;  // threshold for stopping the probabilistic algorithm
                                                 // for computing state accessibility
  const int ACCESSIBILITY_LENGTH = 100;  // maximum sequence length for the probabilistic algorithm
                                         // for computing state accessibility

  const double NOISE_PROBABILITY = 0.05;  // perturbation of observation probabilities
  const double MEAN_SHIFT_COEFF = 0.1;   // coefficient for shifting continuous observation distributions

  const int MIN_NB_ELEMENT = 10;         // minimum size of the sample built by rounding
  const int OBSERVATION_COEFF = 10;      // rounding coefficient for the parametric observation distribution estimator

  const int GAMMA_MAX_NB_DECIMAL = 6;     // maximum number of decimals for the simulation of a gamma distribution
  const int INVERSE_GAUSSIAN_MAX_NB_DECIMAL = 6;  // maximum number of decimals for the simulation
                                                  // of an inverse Gaussian distribution
  const int GAUSSIAN_MAX_NB_DECIMAL = 6;  // maximum number of decimals for the simulation of a Gaussian distribution
  const int DEGREE_DECIMAL_SCALE = 10;   // factor for determining the number of decimals
                                         // for the simulation of a von Mises distribution in degrees
  const int RADIAN_DECIMAL_SCALE = 1000;  // factor for determining the number of decimals
                                          // for the simulation of a von Mises distribution in radians

  // const double SELF_TRANSITION = 0.9;    initial self-tranistion

  enum state_type {
    TRANSIENT ,
    RECURRENT ,
    ABSORBING
  };

  enum model_type {
    MIXTURE ,
    HIDDEN_MARKOV
  };

  enum observation_process {
    CATEGORICAL_PROCESS ,
    DISCRETE_PARAMETRIC ,
    CONTINUOUS_PARAMETRIC ,
    DEFAULT_PROCESS
  };

  enum count_pattern {
    RUN ,
    OCCURRENCE
  };

  enum latent_structure_algorithm {
    NO_LATENT_STRUCTURE ,
    FORWARD ,
//    FORWARD_BACKWARD ,
    VITERBI ,
//    VITERBI_FORWARD_BACKWARD ,
    GENERALIZED_VITERBI ,
    FORWARD_BACKWARD_SAMPLING ,
//    GIBBS_SAMPLING ,
    FORWARD_DYNAMIC_PROGRAMMING
  };



/****************************************************************
 *
 *  Class definition
 */


  class ChainData;

  /// \brief Markov chain

  class Chain {

    friend std::ostream& operator<<(std::ostream &os , const Chain &chain)
    { return chain.ascii_print(os , true); }

  public :

    process_type type;      ///< process type (ORDINARY/EQUILIBRIUM)
    int nb_state;           ///< number of states
    int nb_row;             ///< number of rows of the transition probability matrix
    bool **accessibility;   ///< state accessibility matrix
    int nb_component;       ///< number of classes
    int *component_nb_state;  ///< numbers of states per class
    int **component;        ///< classes
    state_type *stype;      ///< state types (TRANSIENT/RECURRENT/ABSORBING)
    double *initial;        ///< initial probabilities
    double *cumul_initial;  ///< cumulative initial distribution function
    double **transition;    ///< transition probability matrix
    double **cumul_transition;  ///< cumulative transition distribution functions

    void parameter_copy(const Chain&);
    void copy(const Chain&);
    void remove();

    Chain(process_type itype = ORDINARY , int inb_state = 0 , bool init_flag = true);
    Chain(process_type itype , int inb_state , int inb_row , bool init_flag);
    Chain(const Chain &chain) { copy(chain); }
    ~Chain();
    Chain& operator=(const Chain &chain);

    static Chain* parsing(StatError &error , ifstream &in_file , int &line , process_type type);

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;

    void create_cumul();
    void cumul_computation();
    void remove_cumul();
    void log_computation();

    bool** logic_transition_computation() const;
    bool strongly_connected_component_research(StatError &error , bool **ilogic_transition = NULL) const;
    void graph_accessibility_computation(bool **ilogic_transition);
    void probability_accessibility_computation();
    void component_computation(bool **ilogic_transition = NULL);
    bool parallel_initial_state() const;

    void thresholding(double min_probability , bool semi_markov = false);

    int nb_parameter_computation(double min_probability = 0.) const;
    double chi2_value_computation(const ChainData &chain_data) const;

    void init(bool left_right , double self_transition);

    double likelihood_computation(const ChainData &chain_data , bool initial_flag = true) const;
    void chi2_fit(const ChainData &chain_data , Test &test) const;
  };


  /// \brief Data structure corresponding to a Markov chain

  class ChainData : public ChainReestimation<int> {
  public :

    ChainData(process_type itype , int inb_state , int inb_row , bool init_flag = false)
    :ChainReestimation<int>(itype , inb_state , inb_row , init_flag) {}
    ChainData(const ChainData &chain_data);

    int nb_parameter_computation() const;

    void estimation(Chain &chain) const;
  };


  /// \brief Categorical observation process

  class CategoricalProcess {

  public :

    int nb_state;           ///< number of states
    int nb_value;           ///< number of categories
    Distribution **observation;  ///< categorical observation distributions
    Distribution *weight;   ///< theoretical weights of observation distributions
    Distribution *mixture;  ///< mixture of observation distributions
    Distribution *restoration_weight;  ///< weights of observation distributions
                                       ///< deduced from the restoration
    Distribution *restoration_mixture;  ///< mixture of observation distributions

    void copy(const CategoricalProcess &process);
    void remove();

    CategoricalProcess(int inb_state = 0 , int inb_value = 0 , bool observation_flag = false);
    CategoricalProcess(int inb_state , int inb_value , double **observation_probability);
    CategoricalProcess(int inb_state , Distribution **pobservation);
    CategoricalProcess(const CategoricalProcess &process)
    { copy(process); }
    ~CategoricalProcess();
    CategoricalProcess& operator=(const CategoricalProcess &process);

    static CategoricalProcess* parsing(StatError &error , ifstream &in_file ,
                                       int &line , int nb_state ,
                                       model_type model , bool hidden);
    static CategoricalProcess** old_parsing(StatError &error , ifstream &in_file ,
                                            int &line , int nb_state , int &nb_output_process);

    std::ostream& ascii_print(std::ostream &os , FrequencyDistribution **empirical_observation ,
                              FrequencyDistribution *marginal_distribution ,
                              bool exhaustive , bool file_flag , model_type model = HIDDEN_MARKOV) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    FrequencyDistribution **empirical_observation = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL ,
                                    model_type model = HIDDEN_MARKOV) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    FrequencyDistribution **empirical_observation = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL ,
                    model_type model = HIDDEN_MARKOV) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        FrequencyDistribution **empirical_observation = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL ,
                        model_type model = HIDDEN_MARKOV) const;

    bool test_hidden() const;
    void thresholding(double min_probability);
    void state_permutation(int *permut) const;  // permutation of states - to be rework with J.-B.
    int nb_parameter_computation(double min_probability) const;
    Distribution* mixture_computation(Distribution *pweight);
    void init();
  };


  /// \brief Discrete parametric observation process

  class DiscreteParametricProcess {

  public :

    int nb_state;           ///< number of states
    int nb_value;           ///< number of values
    DiscreteParametric **observation;  ///< discrete parametric observation distributions
    Distribution *weight;   ///< theoretical weights of observation distributions
    Distribution *mixture;  ///< mixture of observation distributions
    Distribution *restoration_weight;  ///< weights of observation distributions
                                       ///< deduced from the restoration
    Distribution *restoration_mixture;  ///< mixture of observation distributions

    void copy(const DiscreteParametricProcess &process);
    void remove();

    DiscreteParametricProcess(int inb_state = 0 , int inb_value = 0);
    DiscreteParametricProcess(int inb_state , DiscreteParametric **pobservation);
    DiscreteParametricProcess(const DiscreteParametricProcess &process)
    { copy(process); }
    ~DiscreteParametricProcess();
    DiscreteParametricProcess& operator=(const DiscreteParametricProcess &process);

    static DiscreteParametricProcess* parsing(StatError &error , ifstream &in_file ,
                                              int &line , int nb_state , model_type model ,
                                              double cumul_threshold = OBSERVATION_THRESHOLD);

    std::ostream& ascii_print(std::ostream &os , FrequencyDistribution **empirical_observation ,
                              FrequencyDistribution *marginal_distribution ,
                              bool exhaustive , bool file_flag , model_type model = HIDDEN_MARKOV) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    FrequencyDistribution **empirical_observation = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL ,
                                    model_type model = HIDDEN_MARKOV) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    FrequencyDistribution **empirical_observation = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL ,
                    model_type model = HIDDEN_MARKOV) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        FrequencyDistribution **empirical_observation = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL ,
                        model_type model = HIDDEN_MARKOV) const;

    void nb_value_computation();
    void state_permutation(int *permut) const;  // permutation of states - to be rework with J.-B.
    int nb_parameter_computation() const;
    double mean_computation(Distribution *pweight) const;
    double variance_computation(Distribution *pweight , double mean = D_INF) const;
    Distribution* mixture_computation(Distribution *pweight);
    void init();
  };


  /// \brief Continuous parametric observation process

  class ContinuousParametricProcess {

  public :

    int nb_state;           ///< number of states
    continuous_parametric ident;  ///< identifiers of observation distributions
    bool tied_location;     ///< flag tied means  (gamma, Gaussian)
    bool tied_dispersion;   ///< flag tied dispersion parameters (gamma, Gaussian, von Mises)
    double offset;          ///< offset for Gaussian mixture with evenly spaced means
    angle_unit unit;        ///< unit (degree/radian) for von Mises distributions
    ContinuousParametric **observation;  ///< continuous observation distributions
    Distribution *weight;   ///< theoretical weights of observation distributions
    Distribution *restoration_weight;  ///< weights of observation distributions
                                       ///< deduced from the restoration

    void copy(const ContinuousParametricProcess &process);
    void remove();

    ContinuousParametricProcess(int inb_state = 0);
    ContinuousParametricProcess(int inb_state , ContinuousParametric **pobservation);
    ContinuousParametricProcess(const ContinuousParametricProcess &process)
    { copy(process); }
    ~ContinuousParametricProcess();
    ContinuousParametricProcess& operator=(const ContinuousParametricProcess &process);

    static ContinuousParametricProcess* parsing(StatError &error , ifstream &in_file ,
                                                int &line , int nb_state , model_type model ,
                                                continuous_parametric last_ident = VON_MISES);

    std::ostream& ascii_print(std::ostream &os , Histogram **observation_histogram ,
                              FrequencyDistribution **observation_distribution ,
                              Histogram *marginal_histogram ,
                              FrequencyDistribution *marginal_distribution ,
                              bool exhaustive , bool file_flag , model_type model = HIDDEN_MARKOV) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    Histogram **observation_histogram = NULL ,
                                    FrequencyDistribution **observation_distribution = NULL ,
                                    Histogram *marginal_histogram = NULL ,
                                    FrequencyDistribution *marginal_distribution = NULL ,
                                    model_type model = HIDDEN_MARKOV) const;
    bool plot_print(const char *prefix , const char *title ,
                    int process , Histogram **observation_histogram = NULL ,
                    FrequencyDistribution **observation_distribution = NULL ,
                    Histogram *marginal_histogram = NULL ,
                    FrequencyDistribution *marginal_distribution = NULL ,
                    int nb_value = I_DEFAULT , double **empirical_cdf = NULL ,
                    model_type model = HIDDEN_MARKOV) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        Histogram **observation_histogram = NULL ,
                        FrequencyDistribution **observation_distribution = NULL ,
                        Histogram *marginal_histogram = NULL ,
                        FrequencyDistribution *marginal_distribution = NULL ,
                        int nb_value = I_DEFAULT , double **empirical_cdf = NULL ,
                        model_type model = HIDDEN_MARKOV) const;

    int nb_parameter_computation() const;
    double mean_computation(Distribution *pweight) const;
    double variance_computation(Distribution *pweight , double mean = D_INF) const;
    void select_unit(angle_unit iunit);
    void init(continuous_parametric iident , double min_value , double max_value ,
              double mean , double variance);

    std::ostream& interval_computation(std::ostream &os);
  };


  void log_computation(int nb_value , const double *pmass , double *plog);
  double von_mises_concentration_computation(double mean_direction);


};  // namespace stat_tool



# endif
