/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: hidden_markov_tree.h 2722 2007-02-16 14:17:56Z jbdurand $
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
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef HIDDEN_MARKOV_TREE_H
#define HIDDEN_MARKOV_TREE_H

class Distribution;
class DiscreteParametricModel;
class DiscreteDistributionData;

namespace Stat_trees
{

/*! \file hidden_markov_tree.h
    \brief Purpose:
     provide the class "CategoricalTreeProcess" for the representation
        of tree processes controlled by a discrete tree process. Includes
        the characteristic quantity distributions for each random variable and
     the class "HiddenMarkovTree", which corresponds to hidden Markov trees and
     the class "HiddenMarkovTreeData", which combines a "Trees" object
        and a "HiddenMarkovTree"
*/

class HiddenMarkovTreeData;
class HiddenMarkovIndOutTree;

/****************************************************************
 *
 *  Definitions of constants
 */

/// minimal size of the sample built by rounding off
const int TREE_MIN_NB_ELEMENT= 10;
/// rounding off coefficient used for parametric observation
/// distribution estimation
const int TREE_OBSERVATION_COEFF= 10;
/// threshold to stop the EM iterations
const double MARKOV_TREE_LIKELIHOOD_DIFF= 1.e-6;
/// maximal number of iterations for EM
const int MARKOV_TREE_NB_ITER= 500;
/// default number of hidden state trees for the generalized Viterbi algorithm
const int NB_STATE_TREES= 10;

/// different types of entropy profiles
enum {
  UPWARD ,
  DOWNWARD ,
};

/****************************************************************
 *
 *  Class definitions :
 */

/**
   \class CategoricalTreeProcess
   \brief a collection of categorical distributions, mainly used for
          the representation of emission distributions given
          a state variable (e.g. a hidden state)
          and that of the distribution of characteristic quantities
*/

class CategoricalTreeProcess : public CategoricalProcess
{  // a collection of categorical distributions, mainly used for
   // the representation of emission distributions given
   // a state variable (e.g. a hidden state)
   // and that of the distribution of characteristic quantities

   friend class HiddenMarkovTree;
   friend class HiddenMarkovIndOutTree;
   friend class HiddenMarkovTreeData;
   friend class MarkovOutTree;
   friend class MarkovOutTreeData;

   // friend CategoricalTreeProcess* occupancy_parsing(StatError &error , ifstream &in_file ,
   //                                                  int &line , const Chain &chain ,
   //                                                  double cumul_threshold = CUMUL_THRESHOLD);
   friend bool test_hidden(int nb_output_process, CategoricalTreeProcess **process);

private :

   /// distribution of the tree size
   Distribution *size;
   // Distribution *nb_children;            // [distributions of the quantities for
                                            // which the characteristic distributions
                                            // are not invariant - to determine]
   /// distribution of the tree depth
   Distribution *depth;
   /// probability of non occurrence of a given value
   double *no_occurrence;
   /// probability of leaving a given value
   double *leave;
   /// probability of being absorbed by a given value
   double *absorption;
   /// distribution of the path length before the first occurrence
   /// of a given value, starting from the root
   Distribution **first_occurrence_root;
   /// distribution of path length before the first occurrence
   /// of a given value, starting from all the terminal vertices
   Distribution **first_occurrence_leaves;
   /// distribution of the sizes of the homogeneous zones
   /// for a given value of the variable
   Distribution **sojourn_size;
   /// distribution of the number of homogeneous zones in the tree
   /// for a given value of the variable
   Distribution **nb_zones;
   /// distribution of the occurrence number for a given
   /// value of the variable in the whole tree
   Distribution **nb_occurrences;
   // N.B. check whether these Distributions are DiscreteParametric or not

   void create_characteristic(const Distribution& hsize, bool sojourn_time_flag= true,
                              bool counting_flag= true);

   void copy(const CategoricalTreeProcess& process, bool characteristic_flag= true);

   void copy_double_array(double*& dest, const double* source, int inb_value);
   void copy_Distribution_array(Distribution**& dest,
                                Distribution** const source,
                                int inb_value);
   void remove();
   void remove_double_array(double*& d);
   /** Deallocate arrays of pointers on Distributions */
   void remove_Distribution_array(Distribution**& d, int inb_value);

   void init_Distribution_array(Distribution**& d, int inb_value);
   void init_occupancy(const CategoricalTreeProcess& process, int occupancy_nb_value);

   /** Permutation of the states of \e self */
   void state_permutation(int* perm) const;

   double get_double_characteristic(double* d, int value) const;
   Distribution** get_ptptDistribution_characteristic(Distribution** d) const;
   Distribution* get_ptDistribution_characteristic(Distribution** d,
                                                   int value) const;

   std::ostream& ascii_print(std::ostream &os, int process,
                             FrequencyDistribution** empirical_observation,
                             const TreeCharacteristics* characteristics ,
                             bool exhaustive, bool file_flag) const;
   std::ostream& spreadsheet_print(std::ostream &os, int process,
                                   FrequencyDistribution** empirical_observation= NULL,
                                   const TreeCharacteristics* characteristics= NULL) const;

   /** Gnuplot output of \e self*/
   bool plot_print(const char * prefix, const char * title, int process,
                   FrequencyDistribution** empirical_observation= NULL,
                   const TreeCharacteristics * characteristics= NULL,
                   const FrequencyDistribution * hsize= NULL) const;

   /** Matplotlib output of CategoricalTreeProcess */
   MultiPlotSet* plotable_write(MultiPlotSet &plot, int &index,
                                int process, FrequencyDistribution * const * empirical_observation = NULL,
                                const TreeCharacteristics * characteristics = NULL,
                                const FrequencyDistribution * hsize = NULL) const;

   /** Return the number of views (i.e. the size) in Matplotlib output */
   unsigned int nb_plot_set_computation(int process, FrequencyDistribution * const * empirical_observation = NULL,
                                        const TreeCharacteristics * characteristics = NULL,
                                        const FrequencyDistribution * hsize = NULL) const;

public :

   CategoricalTreeProcess(int inb_state= 0, int inb_value= 0,
                          int observation_flag= false);
   // CategoricalTreeProcess(int inb_state , Distribution **occupancy);
   CategoricalTreeProcess(const CategoricalProcess& process);
   CategoricalTreeProcess(const CategoricalTreeProcess& process ,
                          char manip= 'c', int param= true);
   ~CategoricalTreeProcess();
   CategoricalTreeProcess& operator=(const CategoricalTreeProcess& process);

   // access to class members

   Distribution* get_size() const;
   // Distribution* get_nb_children() const; etc.  // [quantities for which the characteristic
                                                   // distributions are invariant - if any -
                                                   // to determine]
   double get_no_occurrence(int value) const;
   double get_leave(int value) const;
   double get_absorption(int value) const;

   /** Permutation of the state of \e self*/
   void state_permutation(StatError& error, int* perm) const;

   Distribution** get_first_occurrence_root() const;
   Distribution* get_first_occurrence_root(int value) const;
   Distribution** get_first_occurrence_leaves() const;
   Distribution* get_first_occurrence_leaves(int value) const;
   Distribution** get_sojourn_size() const;
   Distribution* get_sojourn_size(int value) const;
   Distribution** get_nb_zones() const;
   Distribution* get_nb_zones(int value) const;
   Distribution** get_nb_occurrences() const;
   Distribution* get_nb_occurrences(int value) const;
};

bool test_hidden(int nb_output_process, CategoricalTreeProcess **process);

/**
   \class HiddenMarkovTree
   \brief the common part of all hidden Markov trees
*/
class HiddenMarkovTree : public StatInterface , protected Chain
{  // hidden Markov tree

   friend class HiddenMarkovTreeData;
   friend class HiddenMarkovIndOutTree;

   friend HiddenMarkovTree* hidden_markov_tree_ascii_read(StatError& error,
                                                          const char * path,
                                                          int size, bool counting_flag,
                                                          double cumul_threshold);

   friend std::ostream& operator<<(std::ostream &os, const HiddenMarkovTree& markov);

public :

   typedef double** double_array_2d;
   typedef double*** double_array_3d;
   typedef double**** double_array_4d;

protected :

   /// data associated with the model
   HiddenMarkovTreeData *markov_data;
   /// (maximal) number of children involved in the transition
   /// probability matrix
   int _ch_order;
   /// probabilities to stay in a given state
   int *self_row;
   /// number of integer observation processes
   int _nb_ioutput_process;
   /// number of double (parametric) observation processes
   int _nb_doutput_process;
   /// non-parametric observation processes
   /// npprocess[0]: (hidden) state process
   CategoricalTreeProcess **npprocess;
   /// integer parametric observation processes
   /// pirocess[0]: unused (NULL)
   DiscreteParametricProcess **piprocess;
   /// double (i.e. floating) parametric observation processes
   /// pdrocess[0]: unused (NULL)
   DiscreteParametricProcess **pdprocess;

   /* HiddenMarkovTree(const Chain * pchain, int inb_nonparam, int inb_param,
                     CategoricalProcess** pobservation, int size);
   HiddenMarkovTree(const HiddenMarkovTree& markov,
                    char manip,
                    int param);
   */

   void copy(const HiddenMarkovTree& markov, bool data_flag= true,
             bool characteristic_flag= true);
   void remove();

   virtual std::ostream& ascii_write(std::ostream& os,
                                     const HiddenMarkovTreeData* otrees,
                                     bool exhaustive= false, bool file_flag= false,
                                     const Test * test= NULL, bool ch_order_flag= true) const;
   virtual std::ostream& spreadsheet_write(std::ostream& os, const HiddenMarkovTreeData * tree,
                                           const Test * test= NULL) const;
   virtual bool plot_write(const char * prefix, const char * title,
                           const HiddenMarkovTreeData * otrees) const;

   void self_row_computation();  // ?
   void component_computation(); // ?

   // void index_state_distribution();
   double* memory_computation() const;
   bool state_profile_write(StatError& error, std::ostream& os,
                            const HiddenMarkovTreeData& tree,
                            int identifier= I_DEFAULT, char format= 'a',
                            int algorithm= UPWARD) const;
   /** Write Gnuplot output of state and entropy profiles */
   bool state_profile_plot_write(StatError& error, const char* prefix ,
                                 const HiddenMarkovTreeData& tree, int identifier,
                                 int vertex, const char* title= NULL,
                                 int algorithm= UPWARD) const;

   void init(bool left_right, double self_transition);
   void init(double self_transition);
   void log_computation();

   // virtual functions common to all hidden_markov_trees ?
   virtual void state_no_occurrence_probability(int state, double increment= LEAVE_INCREMENT);
   virtual void state_first_occurrence_root_distribution(int state, int min_nb_value= 1,
                                                 double cumul_threshold= CUMUL_THRESHOLD);
   virtual void state_first_occurrence_leaves_distribution(int state, int min_nb_value= 1,
                                                   double cumul_threshold= CUMUL_THRESHOLD);
   virtual void state_leave_probability(const double * memory, int state,
                                        double increment= LEAVE_INCREMENT);
   virtual void state_sojourn_size_distribution(const double * memory, int state,
                                                int min_nb_value= 1,
                                                double cumul_threshold= CUMUL_THRESHOLD);
   virtual void state_nb_pattern_mixture(int state, char pattern);

   virtual void output_no_occurrence_probability(int variable, int output,
                                         double increment= LEAVE_INCREMENT);
   virtual void output_first_occurrence_root_distribution(int variable, int output,
                                                  int min_nb_value= 1,
                                                  double cumul_threshold= CUMUL_THRESHOLD);
   virtual void output_first_occurrence_leaves_distribution(int variable, int output,
                                                    int min_nb_value= 1,
                                                    double cumul_threshold= CUMUL_THRESHOLD);
   virtual void output_leave_probability(const double * memory,
                                 int variable, int output,
                                 double increment= LEAVE_INCREMENT);
   virtual void output_sojourn_size_distribution(const double * memory, int variable,
                                         int output, int min_nb_value= 1 ,
                                         double cumul_threshold= CUMUL_THRESHOLD);
   virtual void output_nb_zones_mixture(int variable, int output);
   virtual void output_nb_occurrences_mixture(int variable, int output);

   /**Compute the number of independent parameters of \e self */
   virtual int nb_parameter_computation(double min_probability= 0.) const;
   /**Compute an adaptive penalty for \e self used for model selection*/
   double penalty_computation(double min_probability= 0.) const;

   virtual double_array_3d state_marginal_distribution(const HiddenMarkovTreeData& trees) const;
   virtual double** state_marginal_distribution(const Trees& trees,
                                                int index) const;

   virtual double likelihood_correction(const HiddenMarkovTreeData& trees) const;

   /**Compute the conditional distributions of the observations*/
   void output_conditional_distribution(const HiddenMarkovTreeData& trees,
                                        double_array_3d& output_cond,
                                        bool logcomputation= false,
                                        int index= I_DEFAULT) const;

   /**Compute the maximal marginal entropy, the state tree entropy
   and return the completed likelihood */
   virtual double upward_downward(const HiddenMarkovTreeData& trees,
                                  double& max_marginal_entropy,
                                  double& entropy1,
                                  double& likelihood,
                                  std::deque<int>*& vd,
                                  int index= I_DEFAULT, std::ostream* os= NULL,
                                  char format= 'a',
                                  int vertex= I_DEFAULT,
                                  int entropy_algo=UPWARD) const;

   /** Compute the entropy of partial state processes */
   double partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                      int t,
                                      double_array_3d output_cond_prob,
                                      double_array_3d marginal_prob,
                                      double_array_3d upward_parent_prob,
                                      double_array_3d downward_prob,
                                      double_array_4d downward_pair_prob,
                                      double_array_3d state_entropy,
                                      double_array_3d conditional_entropy,
                                      double_array_4d conditional_prob,
                                      double_array_2d& partial_entropy,
                                      int entropy_algo= UPWARD) const;

   /** Compute the entropy of partial state processes by a downward algorithm*/
   virtual double downward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                                       int t, double_array_3d output_cond_prob,
                                                       double_array_3d marginal_prob,
                                                       double_array_3d upward_parent_prob,
                                                       double_array_3d downward_prob,
                                                       double_array_4d downward_pair_prob,
                                                       double_array_3d state_entropy,
                                                       double_array_3d conditional_entropy,
                                                       double_array_4d conditional_prob,
                                                       double_array_2d& partial_entropy) const;

   /** Compute the entropy of partial state processes by an upward algorithm*/
   virtual double upward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                                     int t,
                                                     double_array_3d downward_prob,
                                                     double_array_3d state_entropy,
                                                     double*& partial_entropy) const;

   /** Compute the conditional entropy*/
   void conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                        double_array_3d marginal_prob,
                                        double_array_3d upward_prob,
                                        double_array_3d upward_parent_prob,
                                        double_array_3d downward_prob,
                                        double_array_2d& expected_conditional_entropy,
                                        double_array_3d& conditional_entropy,
                                        double_array_4d& conditional_prob,
                                        double_array_3d& state_entropy,
                                        int index= I_DEFAULT,
                                        int entropy_algo= UPWARD) const;

   /** Compute the conditional entropy by an upward algorithm*/
   virtual double upward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                         double_array_3d marginal_prob,
                                                         double_array_3d upward_prob,
                                                         double_array_3d upward_parent_prob,
                                                         double_array_3d downward_prob,
                                                         double_array_2d& conditional_entropy,
                                                         int index= I_DEFAULT) const;

   /** Compute the conditional entropy by a downward algorithm*/
   virtual double downward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                           double_array_3d marginal_prob,
                                                           double_array_3d downward_prob,
                                                           double_array_3d upward_prob,
                                                           double_array_3d upward_parent_prob,
                                                           double_array_2d& expected_conditional_entropy,
                                                           double_array_3d& conditional_entropy,
                                                           double_array_4d& conditional_prob,
                                                           double_array_3d& state_entropy,
                                                           int index= I_DEFAULT) const;

   /**Compute the smoothed probabilities and return the likelihood */
   virtual double smoothed_probabilities(const HiddenMarkovTreeData& trees,
                                         double_array_3d& smoothed_prob,
                                         double_array_2d& marginal_entropy,
                                         double_array_2d& conditional_entropy,
                                         double_array_2d& partial_entropy,
                                         int index,
                                         int entropy_algo= UPWARD) const;

   /**Compute the maximal completed likelihood */
   virtual double viterbi(const HiddenMarkovTreeData& trees,
                          int index= I_DEFAULT) const;

   /**Compute the number of possible state trees */
   virtual long double nb_state_trees(const HiddenMarkovTreeData& trees,
                                      int index= I_DEFAULT) const;

   /**Upward_downward Viterbi algorithm */
   virtual HiddenMarkovTreeData* viterbi_upward_downward(const HiddenMarkovTreeData& trees,
                                                         std::vector<ostringstream*>& messages,
                                                         double likelihood,
                                                         double& state_likelihood,
                                                         std::deque<int>*& vd,
                                                         int index= I_DEFAULT,
                                                         std::ostream* os= NULL,
                                                         char format= 'a',
                                                         int vertex= I_DEFAULT) const;
   /**Upward_downward Viterbi algorithm */
   virtual HiddenMarkovTreeData* generalized_viterbi(const HiddenMarkovTreeData& trees,
                                                     std::vector<ostringstream*>& messages,
                                                     int nb_state_trees,
                                                     double likelihood,
                                                     int index,
                                                     int sroot= I_DEFAULT) const;

   /**Compute the completed likelihood for a given set of trees and state trees*/
   virtual double state_likelihood_computation(const HiddenMarkovTreeData& trees) const;

   /**Compute the completed likelihood for one tree and one state tree*/
   virtual double state_likelihood_computation(const HiddenMarkovTreeData& trees,
                                               int index) const;

public :

   HiddenMarkovTree();
   /** Create a Hidden Markov Tree from the type of transition probability matrix,
   the number of states, the number of integer and floating observation processes,
   the number of values for each integer process
   and an indicator array for forcing parametric integer processes */
   HiddenMarkovTree(char itype, int inb_state, int ich_order,
                    int inb_ioutput_process, int inb_doutput_process,
                    int* nb_value, bool* force_param= NULL);
   HiddenMarkovTree(const Chain * pchain, int ich_order, int inb_ioutput_process,
                    CategoricalProcess** pobservation,
                    int size, bool counting_flag);
   HiddenMarkovTree(const Chain * pchain,
                    int ich_order, int inb_ioutput_process, int inb_doutput_process,
                    CategoricalProcess** categorical_observation,
                    DiscreteParametricProcess** iparametric_observation,
                    DiscreteParametricProcess** dparametric_observation,
                    int size, bool counting_flag);
   HiddenMarkovTree(const HiddenMarkovTree& markov, bool data_flag= true,
                    bool characteristic_flag= true);
   virtual ~HiddenMarkovTree();
   HiddenMarkovTree& operator=(const HiddenMarkovTree& markov);

   /** Return a copy of \e self with the same dynamic class */
   virtual HiddenMarkovTree* HiddenMarkovTreeCopy(bool data_flag = true,
                                                  bool characteristic_flag = true) const;

   DiscreteParametricModel* extract(StatError& error, int type,
                                    int variable, int value) const;

   /** Return the data part of a HiddenMarkovTree,
       keeping a reference on \e self */
   virtual HiddenMarkovTreeData* extract_data(StatError& error) const;

   HiddenMarkovTree* thresholding(double min_probability= MIN_PROBABILITY) const;

   /** Permutation of the states of \e self */
   void state_permutation(StatError& error, int* perm) const;

   virtual std::ostream& line_write(std::ostream& os) const;

   virtual std::ostream& ascii_write(std::ostream& os, bool exhaustive= false) const;
   virtual bool ascii_write(StatError& error, const char * path,
                            bool exhaustive= false) const;
   virtual bool spreadsheet_write(StatError& error, const char * path) const;
   virtual bool plot_write(StatError& error, const char * prefix,
                           const char * title= NULL) const;

   void characteristic_computation(int size, bool counting_flag,
                                   int variable= I_DEFAULT);
   void characteristic_computation(const HiddenMarkovTreeData& tree,
                                   bool counting_flag,
                                   int variable= I_DEFAULT,
                                   bool size_flag= true);

   /** Compute the state profiles for a given tree,
       including the smoothed probabilities,
       the upward-downward and the generalized Viterbi algorithms */
   bool state_profile(StatError& error,
                      const HiddenMarkovTreeData& tree,
                      int index,
                      HiddenMarkovTreeData*& smoothed_state_tree,
                      HiddenMarkovTreeData*& nstate_trees,
                      HiddenMarkovTreeData*& viterbi_upward_downward,
                      HiddenMarkovTreeData*& generalized_restoration,
                      std::vector<ostringstream*>& messages,
                      int state_tree= GENERALIZED_VITERBI,
                      unsigned int nb_state_trees= NB_STATE_TREES,
                      int entropy_algo= UPWARD,
                      int root= I_DEFAULT) const;

   /** Write Gnuplot output of state and entropy profiles */
   virtual bool tree_state_profile_plot_write(StatError &error, const char *prefix,
                                              const HiddenMarkovTreeData& trees,
                                              int identifier, int vertex,
                                              const char *title= NULL,
                                              int entropy_algo= UPWARD) const;

   /** Write Gnuplot output of state and entropy profiles */
   bool state_profile_plot_write(StatError &error, const char *prefix,
                                 int identifier, int vertex,
                                 const char *title= NULL,
                                 int entropy_algo= UPWARD) const;


   // virtual functions common to all hidden_markov_trees ?

   virtual double likelihood_computation(const Trees& trees, int index) const;
   virtual double likelihood_computation(HiddenMarkovTreeData& trees) const;

   /** Compute optimal state trees */
   virtual HiddenMarkovTreeData* state_tree_computation(StatError& error,
                                                        const Trees& trees,
                                                        int algorithm= VITERBI,
                                                        bool characteristic_flag= false,
                                                        int index= I_DEFAULT) const;

    // access to class members

   /** Return the data part of a HiddenMarkovTree,
       losing a reference on \e self */
   HiddenMarkovTreeData* get_markov_data() const;
   int get_ch_order() const;
   /** return the number of states of \e self*/
   int get_nb_state() const;
   int get_self_row(int state) const;
   /** return the number of processes with discrete values*/
   int get_nb_ioutput_process() const;
   /** return the number of floating processes*/
   int get_nb_doutput_process() const;
   /** return the number of values of a process with finite values*/
   int get_nb_values(int variable) const;
   /** return "true" if process ivariable is parametric */
   bool is_parametric(int ivariable) const;
   CategoricalTreeProcess** get_categorical_process() const;
   CategoricalTreeProcess* get_categorical_process(int variable) const;
   DiscreteParametricProcess** get_iparametric_process() const;
   DiscreteParametricProcess* get_iparametric_process(int variable) const;
   DiscreteParametricProcess** get_dparametric_process() const;
   DiscreteParametricProcess* get_dparametric_process(int variable) const;
};

HiddenMarkovTree* hidden_markov_tree_ascii_read(StatError& error,
                                                const char * path,
                                                int size= I_DEFAULT_TREE_SIZE,
                                                bool counting_flag= true,
                                                double cumul_threshold= OCCUPANCY_THRESHOLD);
/**
   \class HiddenMarkovTreeData
   \brief set of tree-structured data potentially associated with
          a statistical model having an underlying structure (states)
*/

class HiddenMarkovTreeData : public Trees
{  // a pair (HiddenMarkovTree model, Trees)

   // friend classes
   friend class HiddenMarkovTree;
   friend class HiddenMarkovIndOutTree;
   // Above lines may require to be added to Typed_edge_trees

   friend std::ostream& operator<<(std::ostream &os, const HiddenMarkovTreeData& trees);

public :

   typedef Trees::tree_type tree_type;
   typedef tree_type::value value;

   typedef Typed_edge_one_int_tree::tree_type state_tree_type;
   typedef state_tree_type::value state_value;

   typedef tree_type::key key;
   typedef tree_type::vertices_size_type vertices_size_type;
   typedef tree_type::children_iterator children_iterator;
   typedef tree_type::vertex_iterator vertex_iterator;

   typedef FrequencyDistribution*** ptFrequencyDistribution_array_2d;
   typedef TreeCharacteristics::ptFrequencyDistribution_array ptFrequencyDistribution_array;
   typedef Histogram*** ptHistogram_array_2d;
   typedef TreeCharacteristics::ptHistogram_array ptHistogram_array;
   typedef Typed_edge_one_int_tree** ptOne_int_tree_set;
   typedef std::vector<HiddenMarkovTreeData*> pt_hmtd_vector;
   typedef double** double_array_2d;

private :

   /// the hidden Markov tree model
   HiddenMarkovTree* markov;
   /// initial state and transitions
   ChainData* chain_data;
   /// warning: the names of the following two private members
   /// are switched with respect to STAT::Hidden_markov_data
   /// likelihood for the given trees
   double likelihood;
   /// completed likelihood
   double hidden_likelihood;
   /// number of hidden states
   int _nb_states;
   /// entropy of state trees given observations
   double sample_entropy;
   /// entropy of each state tree given observed tree
   double* entropy;

   /// hidden trees
   ptOne_int_tree_set state_trees;

   /// FrequencyDistribution corresponding to the conditional observation distribution,
   /// depending on the considered observed (integral) variable
   /// and on the value of the state variable
   ptFrequencyDistribution_array_2d observation_distribution;

   /// Histogram corresponding to the conditional observation distribution,
   /// depending on the considered observed (continuous) variable
   /// and on the value of the state variable
   ptHistogram_array_2d observation_histogram;

    /// frequency distribution of the characteristic quantities
    /// for the hidden state variable
   TreeCharacteristics *state_characteristics;

   void copy(const HiddenMarkovTreeData& otrees, bool model_flag= true,
             bool characteristic_flag= true);

   ///  Deallocation. The markov part is deleted if any.
   void remove();

   void nb_state_computation();
   void observation_frequency_distribution_computation(int ivariable);
   // computation of the frequency distributions for the conditional distribution
   // of one observed variable
   void build_state_trees();
   void build_state_characteristics();

   /** Return a HiddenMarkovTreeData without the state variable */
   HiddenMarkovTreeData* remove_state_variable() const;

   /** Permutation of the states of \e self */
   void state_permutation(int* perm);

   std::ostream& state_profile_ascii_print(std::ostream& os, int index,
                                           int nb_state,
                                           double_array_2d smoothed) const;

   std::ostream& state_profile_spreadsheet_print(std::ostream& os, int index,
                                                 int nb_state,
                                                 double_array_2d smoothed) const;

   std::ostream& state_profile_plot_print(std::ostream& os, int index,
                                          int nb_state,
                                          double_array_2d smoothed) const;

   /** Write Gnuplot output of state and entropy profiles */
   std::ostream& profile_plot_print(std::ostream& os, int index, int nb_state,
                                    double_array_2d smoothed,
                                    double_array conditional_entropy,
                                    double_array marginal_entropy,
                                    double_array partial_entropy,
                                    key vertex,
                                    generic_visitor<tree_type>::vertex_deque*& path) const;

   /** Write Gnuplot output of Viterbi state profiles */
   std::ostream& profile_plot_print(std::ostream& os, int index, int nb_state,
                                    double_array_2d ratio,
                                    key vertex,
                                    generic_visitor<tree_type>::vertex_deque*& path) const;

public :

   HiddenMarkovTreeData();
   HiddenMarkovTreeData(int inb_integral,
                        int inb_float,
                        const FrequencyDistribution& ihsize,
                        const FrequencyDistribution& ihnb_children,
                        bool no_child_flag,
                        bool init_flag);
   HiddenMarkovTreeData(int inb_integral, int inb_float, int inb_trees);
   HiddenMarkovTreeData(int inb_trees,
                        int* itype,
                        Default_tree** otrees);
   HiddenMarkovTreeData(const Trees& otrees);
   HiddenMarkovTreeData(const HiddenMarkovTreeData& trees,
                        bool model_flag= true,
                        bool characteristic_flag= true);
   ~HiddenMarkovTreeData();
   HiddenMarkovTreeData& operator=(const HiddenMarkovTreeData& trees);

   DiscreteDistributionData* extract(StatError& error, int type,
                                     int variable, int value) const;

   /** Return the model part of a HiddenMarkovTreeData,
       keeping a reference on \e self */
   HiddenMarkovTree* extract_model(StatError& error) const;

   /** Return mixture distribution with frequency distribution for a given variable */
   DiscreteMixtureData* extract_marginal(StatError& error, int variable) const;

   /*
   Print mixture distribution with frequency distribution for a given variable
   bool marginal_ascii_write(StatError& error, int variable,
                             const char *path, bool exhaustive= true) const; */

   HiddenMarkovTreeData* merge(StatError& error,
                               const pt_hmtd_vector& otrees) const;

   /* Print \e self on a single line
   std::ostream& line_write(std::ostream& os) const; */

   std::ostream& ascii_write(std::ostream& os, bool exhaustive= false) const;
   bool ascii_write(StatError& error, const char *path,
                    bool exhaustive= false) const;
   bool spreadsheet_write(StatError& error, const char * path) const;
   bool plot_write(StatError& error, const char * prefix,
                   const char * title= NULL) const;

   // model identification

   /** Estimate a Hidden Markov Out Tree from an initial model and the data in \e self */
   HiddenMarkovIndOutTree* hidden_markov_ind_out_tree_estimation(StatError& error,
                                                          std::ostream& os,
                                                          const HiddenMarkovIndOutTree& ihmarkov,
                                                          bool counting_flag= true,
                                                          int state_trees= VITERBI,
                                                          int algorithm= FORWARD_BACKWARD,
                                                          double saem_exponent= 1.,
                                                          int nb_iter= I_DEFAULT,
                                                          bool force_param= false) const;

/*   HiddenMarkovIndOutTree* hidden_markov_ind_out_tree_estimation2(StatError& error,
                                                           std::ostream& os,
                                                           const HiddenMarkovIndOutTree& ihmarkov,
                                                           bool counting_flag= true,
                                                           int state_trees= VITERBI,
                                                           int algorithm= FORWARD_BACKWARD,
                                                           double saem_exponent= 1.,
                                                           int nb_iter= I_DEFAULT,
                                                           bool force_param= false) const;*/

   /** Estimate a Hidden Markov Out Tree from the number of states and the data in \e self */
   HiddenMarkovIndOutTree* hidden_markov_ind_out_tree_estimation(StatError& error,
                                                          std::ostream& os,
                                                          char type,
                                                          int nb_state,
                                                          bool left_right,
                                                          bool counting_flag= true,
                                                          int state_trees= VITERBI,
                                                          int algorithm= FORWARD_BACKWARD,
                                                          double saem_exponent= 1.,
                                                          double self_transition= D_DEFAULT,
                                                          int nb_iter= I_DEFAULT,
                                                          bool* force_param= NULL) const;

   void create_observation_frequency_distribution(int nb_state);
   // allocation of the frequency distributions corresponding to the observation distributions
   void observation_frequency_distribution_computation();
   // computation of the frequency distributions corresponding to the observation distributions
   void build_observation_frequency_distribution();
   // allocation and computation of the frequency distributions
   // corresponding to the observation distributions

   // access to class members

   HiddenMarkovTree* get_markov() const; // { return markov; }
   ChainData* get_chain_data() const; // { return chain_data; }
   double get_likelihood() const; // { return likelihood; }
   double get_hidden_likelihood() const; // { return hidden_likelihood; }
   int get_nb_states() const;

   ptFrequencyDistribution_array_2d get_observation() const;
   ptFrequencyDistribution_array get_observation(int variable) const;
   FrequencyDistribution* get_observation(int variable, int state) const;

   TreeCharacteristics* get_state_characteristics() const;
   // access to the characteristic quantity distributions
   // for the state variable

   /** Return the set of state trees (a new instance is allocated) */
   ptOne_int_tree_set get_state_trees() const;
   /** Return the set of state trees (return a pointer; object should not be deallocated) */
   ptOne_int_tree_set get_state_trees_ptr() const;
   /** Return a given state tree (a new instance is allocated) */
   Typed_edge_one_int_tree* get_state_tree(int itree) const;
   /** Return a given state tree (return a pointer; object should not be deallocated) */
   Typed_edge_one_int_tree* get_state_tree_ptr(int itree) const;

   /** Return a HiddenMarkovTreeData containing the states
       as a variable */
   HiddenMarkovTreeData* get_state_hidden_markov_tree_data() const;
   /** For a given tree, and using a given type of entropy profile,
       return a HiddenMarkovTreeData containing the states,
       the entropy and the smoothed probabilities as variables */
   HiddenMarkovTreeData*
      get_state_smoothed_hidden_markov_tree_data(int index= I_DEFAULT,
                                                 int algorithm= VITERBI,
                                                 int entropy_algo= UPWARD) const;
   // access to the state trees

   std::ostream& ascii_write_observation(std::ostream &os,
                                         bool exhaustive,
                                         bool file_flag) const;

};
};

#endif
