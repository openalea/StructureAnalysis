/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/STAT_TREES/src/hidden_markov_tree.h,v $
 *       $Id: hidden_markov_tree.h 2722 2007-02-16 14:17:56Z jbdurand $
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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
class Parametric_model;
class Distribution_data;

namespace Stat_trees
{

/*! \file hidden_markov_tree.h
    \brief Purpose:
     provide the class "Non_parametric_tree_process" for the representation
        of tree processes controlled by a discrete tree process. Includes
        the characteristic quantity distributions for each random variable and
     the class "Hidden_markov_tree", which corresponds to hidden Markov trees and
     the class "Hidden_markov_tree_data", which combines a "Trees" object
        and a "Hidden_markov_tree"
*/

class Hidden_markov_tree_data;

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
   \class Nonparametric_tree_process
   \brief a collection of nonparametric distributions, mainly used for
          the representation of emission distributions given
          a state variable (e.g. a hidden state)
          and that of the distribution of characteristic quantities
*/

class Nonparametric_tree_process : public Nonparametric_process
{  // a collection of nonparametric distributions, mainly used for
   // the representation of emission distributions given
   // a state variable (e.g. a hidden state)
   // and that of the distribution of characteristic quantities

   friend class Hidden_markov_tree;
   friend class Hidden_markov_out_tree;
   friend class Hidden_markov_tree_data;

   // friend Nonparametric_tree_process* occupancy_parsing(Format_error &error , ifstream &in_file ,
   //                                                          int &line , const Chain &chain ,
   //                                                          double cumul_threshold = CUMUL_THRESHOLD);
   friend bool test_hidden(int nb_output_process, Nonparametric_tree_process **process);

private :

   /// distribution of the tree size
   Distribution *size;
   // Distribution *nb_children;            // [distributions of the quantities for
                                            // which the characteristic distributions
                                            // are invariant - if any - to determine]
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
   // N.B. check whether these Distributions are Parametric or not

   void create_characteristic(const Distribution& hsize, bool sojourn_time_flag= true,
                              bool counting_flag= true);

   void copy(const Nonparametric_tree_process& process, bool characteristic_flag= true);

   void copy_double_array(double*& dest, const double* source, int inb_value);
   void copy_Distribution_array(Distribution**& dest,
                                Distribution** const source,
                                int inb_value);
   void remove();
   void remove_double_array(double*& d);
   /** Deallocate arrays of pointers on Distributions */
   void remove_Distribution_array(Distribution**& d, int inb_value);

   void init_Distribution_array(Distribution**& d, int inb_value);
   void init_occupancy(const Nonparametric_tree_process& process, int occupancy_nb_value);

   /** Permutation of the states of \e self */
   void state_permutation(int* perm) const;

   double get_double_characteristic(double* d, int value) const;
   Distribution** get_ptptDistribution_characteristic(Distribution** d) const;
   Distribution* get_ptDistribution_characteristic(Distribution** d,
                                                   int value) const;

   std::ostream& ascii_print(std::ostream &os, int process, Histogram** empirical_observation,
                             const Tree_characteristics* characteristics ,
                             bool exhaustive, bool file_flag) const;
   std::ostream& spreadsheet_print(std::ostream &os, int process,
                              Histogram** empirical_observation= NULL,
                              const Tree_characteristics* characteristics= NULL) const;

   /** Gnuplot output of \e self*/
   bool plot_print(const char * prefix, const char * title, int process,
                   Histogram** empirical_observation= NULL,
                   const Tree_characteristics * characteristics= NULL,
                   const Histogram * hsize= NULL) const;

public :

   Nonparametric_tree_process(int inb_state= 0, int inb_value= 0,
                              int observation_flag= false);
   // Nonparametric_tree_process(int inb_state , Distribution **occupancy);
   Nonparametric_tree_process(const Nonparametric_process& process);
   Nonparametric_tree_process(const Nonparametric_tree_process& process ,
                              char manip= 'c', int param= true);
   ~Nonparametric_tree_process();
   Nonparametric_tree_process& operator=(const Nonparametric_tree_process& process);

   // access to class members

   Distribution* get_size() const;
   // Distribution* get_nb_children() const; etc.  // [quantities for which the characteristic
                                                   // distributions are invariant - if any -
                                                   // to determine]
   double get_no_occurrence(int value) const;
   double get_leave(int value) const;
   double get_absorption(int value) const;

   /** Permutation of the state of \e self*/
   void state_permutation(Format_error& error, int* perm) const;

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

bool test_hidden(int nb_output_process, Nonparametric_tree_process **process);

/**
   \class Hidden_markov_tree
   \brief the common part of all hidden Markov trees
*/
class Hidden_markov_tree : public STAT_interface , protected Chain
{  // hidden Markov tree

   friend class Hidden_markov_tree_data;
   friend class Hidden_markov_out_tree;

   friend Hidden_markov_tree* hidden_markov_tree_ascii_read(Format_error& error,
                                                            const char * path,
                                                            int size, bool counting_flag,
                                                            double cumul_threshold);

   friend std::ostream& operator<<(std::ostream &os, const Hidden_markov_tree& markov);

public :

   typedef double** double_array_2d;
   typedef double*** double_array_3d;
   typedef double**** double_array_4d;

protected :

   /// data associated with the model
   Hidden_markov_tree_data *markov_data;
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
   Nonparametric_tree_process **npprocess;
   /// integer parametric observation processes
   Parametric_process **piprocess;
   /// double (i.e. floating) parametric observation processes
   Parametric_process **pdprocess;

   /* Hidden_markov_tree(const Chain * pchain, int inb_nonparam, int inb_param,
                      Nonparametric_process** pobservation, int size);
   Hidden_markov_tree(const Hidden_markov_tree& markov,
                      char manip,
                      int param);
   */

   void copy(const Hidden_markov_tree& markov, bool data_flag= true,
             bool characteristic_flag= true);
   void remove();

   virtual std::ostream& ascii_write(std::ostream& os,
                                     const Hidden_markov_tree_data* otrees,
                                     bool exhaustive= false, bool file_flag= false,
                                     const Test * test= NULL, bool ch_order_flag= true) const;
   virtual std::ostream& spreadsheet_write(std::ostream& os, const Hidden_markov_tree_data * tree,
                                           const Test * test= NULL) const;
   virtual bool plot_write(const char * prefix, const char * title,
                           const Hidden_markov_tree_data * otrees) const;

   void self_row_computation();  // ?
   void component_computation(); // ?

   // void index_state_distribution();
   double* memory_computation() const;
   bool state_profile_write(Format_error& error, std::ostream& os,
                            const Hidden_markov_tree_data& tree,
                            int identifier= I_DEFAULT, char format= 'a',
                            int algorithm= UPWARD) const;
   /** Write Gnuplot output of state and entropy profiles */
   bool state_profile_plot_write(Format_error& error, const char* prefix ,
                                 const Hidden_markov_tree_data& tree, int identifier,
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

   virtual double_array_3d state_marginal_distribution(const Hidden_markov_tree_data& trees) const;
   virtual double** state_marginal_distribution(const Trees& trees,
                                                int index) const;

   virtual double likelihood_correction(const Hidden_markov_tree_data& trees) const;

   /**Compute the conditional distributions of the observations*/
   void output_conditional_distribution(const Hidden_markov_tree_data& trees,
                                        double_array_3d& output_cond,
                                        bool logcomputation= false,
                                        int index= I_DEFAULT) const;

   /**Compute the maximal marginal entropy, the state tree entropy
   and return the completed likelihood */
   virtual double upward_downward(const Hidden_markov_tree_data& trees,
                                  double& max_marginal_entropy,
                                  double& entropy1,
                                  double& likelihood,
                                  std::deque<int>*& vd,
                                  int index= I_DEFAULT, std::ostream* os= NULL,
                                  char format= 'a',
                                  int vertex= I_DEFAULT,
                                  int entropy_algo=UPWARD) const;

   /** Compute the entropy of partial state processes */
   double partial_entropy_computation(const Hidden_markov_tree_data& trees,
                                      int t,
                                      double_array_3d output_cond_prob,
                                      double_array_3d marginal_prob,
                                      double_array_3d upward_parent_prob,
                                      double_array_3d downward_prob,
                                      double_array_3d state_entropy,
                                      double_array_3d conditional_entropy,
                                      double_array_4d conditional_prob,
                                      double*& partial_entropy,
                                      int entropy_algo= UPWARD) const;

   /** Compute the entropy of partial state processes by a downward algorithm*/
   virtual double downward_partial_entropy_computation(const Hidden_markov_tree_data& trees,
                                                       int t,
                                                       double_array_3d output_cond_prob,
                                                       double_array_3d marginal_prob,
                                                       double_array_3d upward_parent_prob,
                                                       double_array_3d downward_prob,
                                                       double_array_3d state_entropy,
                                                       double_array_3d conditional_entropy,
                                                       double_array_4d conditional_prob,
                                                       double*& partial_entropy) const;

   /** Compute the entropy of partial state processes by an upward algorithm*/
   virtual double upward_partial_entropy_computation(const Hidden_markov_tree_data& trees,
                                                     int t,
                                                     double_array_3d downward_prob,
                                                     double_array_3d state_entropy,
                                                     double*& partial_entropy) const;

   /** Compute the conditional entropy*/
   void conditional_entropy_computation(const Hidden_markov_tree_data& trees,
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
   virtual void upward_conditional_entropy_computation(const Hidden_markov_tree_data& trees,
                                                       double_array_3d marginal_prob,
                                                       double_array_3d upward_prob,
                                                       double_array_3d upward_parent_prob,
                                                       double_array_3d downward_prob,
                                                       double_array_2d& conditional_entropy,
                                                       int index= I_DEFAULT) const;

   /** Compute the conditional entropy by a downward algorithm*/
   virtual void downward_conditional_entropy_computation(const Hidden_markov_tree_data& trees,
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
   virtual double smoothed_probabilities(const Hidden_markov_tree_data& trees,
                                         double_array_3d& smoothed_prob,
                                         double_array_2d& marginal_entropy,
                                         double_array_2d& conditional_entropy,
                                         double_array_2d& partial_entropy,
                                         int index,
                                         int entropy_algo= UPWARD) const;

   /**Compute the maximal completed likelihood */
   virtual double viterbi(const Hidden_markov_tree_data& trees,
                          int index= I_DEFAULT) const;

   /**Compute the number of possible state trees */
   virtual long double nb_state_trees(const Hidden_markov_tree_data& trees,
                                      int index= I_DEFAULT) const;

   /**Upward_downward Viterbi algorithm */
   virtual Hidden_markov_tree_data* viterbi_upward_downward(const Hidden_markov_tree_data& trees,
                                                            std::vector<ostringstream*>& messages,
                                                            double likelihood,
                                                            double& state_likelihood,
                                                            std::deque<int>*& vd,
                                                            int index= I_DEFAULT,
                                                            std::ostream* os= NULL,
                                                            char format= 'a',
                                                            int vertex= I_DEFAULT) const;
   /**Upward_downward Viterbi algorithm */
   virtual Hidden_markov_tree_data* generalized_viterbi(const Hidden_markov_tree_data& trees,
                                                        std::vector<ostringstream*>& messages,
                                                        int nb_state_trees,
                                                        double likelihood,
                                                        int index) const;

   /**Compute the completed likelihood for a given set of trees and state trees*/
   virtual double state_likelihood_computation(const Hidden_markov_tree_data& trees) const;

   /**Compute the completed likelihood for one tree and one state tree*/
   virtual double state_likelihood_computation(const Hidden_markov_tree_data& trees,
                                               int index) const;

public :

   Hidden_markov_tree();
   /** Create a Hidden Markov Tree from the type of transition probability matrix,
   the number of states, the number of integer and floating observation processes,
   the number of values for each integer process
   and an indicator array for forcing parametric integer processes */
   Hidden_markov_tree(char itype, int inb_state, int ich_order,
                      int inb_ioutput_process, int inb_doutput_process,
                      int* nb_value, bool* force_param= NULL);
   Hidden_markov_tree(const Chain * pchain, int ich_order, int inb_ioutput_process,
                      Nonparametric_process** pobservation,
                      int size, bool counting_flag);
   Hidden_markov_tree(const Chain * pchain,
                      int ich_order, int inb_ioutput_process, int inb_doutput_process,
                      Nonparametric_process** nonparametric_observation,
                      Parametric_process** iparametric_observation,
                      Parametric_process** dparametric_observation,
                      int size, bool counting_flag);
   Hidden_markov_tree(const Hidden_markov_tree& markov, bool data_flag= true,
                      bool characteristic_flag= true);
   virtual ~Hidden_markov_tree();
   Hidden_markov_tree& operator=(const Hidden_markov_tree& markov);

   Parametric_model* extract(Format_error& error, int type,
                             int variable, int value) const;
   /** Return the data part of a Hidden_markov_tree,
       keeping a reference on \e self */
   virtual Hidden_markov_tree_data* extract_data(Format_error& error) const;

   Hidden_markov_tree* thresholding(double min_probability= MIN_PROBABILITY) const;

   /** Permutation of the states of \e self */
   void state_permutation(Format_error& error, int* perm) const;

   virtual std::ostream& line_write(std::ostream& os) const;

   virtual std::ostream& ascii_write(std::ostream& os, bool exhaustive= false) const;
   virtual bool ascii_write(Format_error& error, const char * path,
                            bool exhaustive= false) const;
   virtual bool spreadsheet_write(Format_error& error, const char * path) const;
   virtual bool plot_write(Format_error& error, const char * prefix,
                           const char * title= NULL) const;

   void characteristic_computation(int size, bool counting_flag,
                                   int variable= I_DEFAULT);
   void characteristic_computation(const Hidden_markov_tree_data& tree,
                                   bool counting_flag,
                                   int variable= I_DEFAULT,
                                   bool size_flag= true);

   /** Compute the state profiles for a given tree,
       including the smoothed probabilities,
       the upward-downward and the generalized Viterbi algorithms */
   bool state_profile(Format_error& error,
                      const Hidden_markov_tree_data& tree,
                      int index,
                      Hidden_markov_tree_data*& smoothed_state_tree,
                      Hidden_markov_tree_data*& nstate_trees,
                      Hidden_markov_tree_data*& viterbi_upward_downward,
                      Hidden_markov_tree_data*& generalized_restoration,
                      std::vector<ostringstream*>& messages,
                      int state_tree= GENERALIZED_VITERBI,
                      unsigned int nb_state_trees= NB_STATE_TREES,
                      int entropy_algo= UPWARD) const;

   /** Write Gnuplot output of state and entropy profiles */
   virtual bool tree_state_profile_plot_write(Format_error &error, const char *prefix,
                                              const Hidden_markov_tree_data& trees,
                                              int identifier, int vertex,
                                              const char *title= NULL,
                                              int entropy_algo= UPWARD) const;

   /** Write Gnuplot output of state and entropy profiles */
   bool state_profile_plot_write(Format_error &error, const char *prefix,
                                 int identifier, int vertex,
                                 const char *title= NULL,
                                 int entropy_algo= UPWARD) const;


   // virtual functions common to all hidden_markov_trees ?

   virtual double likelihood_computation(const Trees& trees, int index) const;
   virtual double likelihood_computation(const Hidden_markov_tree_data& trees) const;

   /** Compute optimal state trees */
   virtual Hidden_markov_tree_data* state_tree_computation(Format_error& error,
                                                           const Trees& trees,
                                                           int algorithm= VITERBI,
                                                           bool characteristic_flag= false,
                                                           int index= I_DEFAULT) const;

    // access to class members

   /** Return the data part of a Hidden_markov_tree,
       losing a reference on \e self */
   Hidden_markov_tree_data* get_markov_data() const;
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
   Nonparametric_tree_process** get_non_parametric_process() const;
   Nonparametric_tree_process* get_non_parametric_process(int variable) const;
   Parametric_process** get_iparametric_process() const;
   Parametric_process* get_iparametric_process(int variable) const;
   Parametric_process** get_dparametric_process() const;
   Parametric_process* get_dparametric_process(int variable) const;
};

Hidden_markov_tree* hidden_markov_tree_ascii_read(Format_error& error,
                                                  const char * path,
                                                  int size= I_DEFAULT_TREE_SIZE,
                                                  bool counting_flag= true,
                                                  double cumul_threshold= OCCUPANCY_THRESHOLD);
/**
   \class Hidden_markov_tree_data
   \brief set of tree-structured data potentially associated with
          a statistical model having an underlying structure (states)
*/

class Hidden_markov_tree_data : public Trees
{  // a pair (Hidden_markov_tree model, Trees)

   // friend classes
   friend class Hidden_markov_tree;
   friend class Hidden_markov_out_tree;
   // Above lines may require to be added to Typed_edge_trees

   friend std::ostream& operator<<(std::ostream &os, const Hidden_markov_tree_data& trees);

public :

   typedef Trees::tree_type tree_type;
   typedef tree_type::value value;

   typedef Typed_edge_one_int_tree::tree_type state_tree_type;
   typedef state_tree_type::value state_value;

   typedef tree_type::key key;
   typedef tree_type::vertices_size_type vertices_size_type;
   typedef tree_type::children_iterator children_iterator;
   typedef tree_type::vertex_iterator vertex_iterator;

   typedef Histogram*** ptHistogram_array_2d;
   typedef Tree_characteristics::ptHistogram_array ptHistogram_array;
   typedef Typed_edge_one_int_tree** ptOne_int_tree_set;
   typedef std::vector<Hidden_markov_tree_data*> pt_hmtd_vector;
   typedef double** double_array_2d;

private :

   /// the hidden Markov tree model
   Hidden_markov_tree* markov;
   /// initial state and transitions
   Chain_data* chain_data;
   /// warning: the names of the following two private members
   /// are switched with respect to STAT::Hidden_markov_data
   /// likelihood for the given trees
   double likelihood;
   /// completed likelihood
   double hidden_likelihood;
   // number of hidden states
   int _nb_states;

   /// hidden trees
   ptOne_int_tree_set state_trees;

   /// histogram corresponding to the conditional observation distribution,
   /// depending on the considered observed (integral) variable
   /// and on the value of the state variable
   ptHistogram_array_2d observation;

    /// histogram of the characteristic quantities
    /// for the hidden state variable
   Tree_characteristics *state_characteristics;

   void copy(const Hidden_markov_tree_data& otrees, bool model_flag= true,
             bool characteristic_flag= true);
   void remove();

   void nb_state_computation();
   void observation_histogram_computation(int ivariable);
   // computation of the histograms for the conditional distribution
   // of one observed variable
   void build_state_trees();
   void build_state_characteristics();

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

   Hidden_markov_tree_data();
   Hidden_markov_tree_data(int inb_integral,
                           int inb_float,
                           const Histogram& ihsize,
                           const Histogram& ihnb_children,
                           bool no_child_flag,
                           bool init_flag);
   Hidden_markov_tree_data(int inb_integral, int inb_float, int inb_trees);
   Hidden_markov_tree_data(int inb_trees,
                           int* itype,
                           Default_tree** otrees);
   Hidden_markov_tree_data(const Trees& otrees);
   Hidden_markov_tree_data(const Hidden_markov_tree_data& trees,
                           bool model_flag= true,
                           bool characteristic_flag= true);
   ~Hidden_markov_tree_data();
   Hidden_markov_tree_data& operator=(const Hidden_markov_tree_data& trees);

   Distribution_data* extract(Format_error& error, int type,
                              int variable, int value) const;

   /*
   Return mixture distribution with histogram for a given variable
   Distribution_data* extract_marginal(Format_error& error, int variable) const;

   Print mixture distribution with histogram for a given variable
   bool marginal_ascii_write(Format_error& error, int variable,
                             const char *path, bool exhaustive= true) const; */

   Hidden_markov_tree_data* merge(Format_error& error,
                                  const pt_hmtd_vector& otrees) const;

   /* Print \e self on a single line
   std::ostream& line_write(std::ostream& os) const; */

   std::ostream& ascii_write(std::ostream& os, bool exhaustive= false) const;
   bool ascii_write(Format_error& error, const char *path,
                    bool exhaustive= false) const;
   bool spreadsheet_write(Format_error& error, const char * path) const;
   bool plot_write(Format_error& error, const char * prefix,
                   const char * title= NULL) const;

   // model identification

   /** Estimate a Hidden Markov Out Tree from an initial model and the data in \e self */
   Hidden_markov_out_tree* hidden_markov_out_tree_estimation(Format_error& error,
                                                             std::ostream& os,
                                                             const Hidden_markov_out_tree& ihmarkov,
                                                             bool counting_flag= true,
                                                             int state_trees= VITERBI,
                                                             int algorithm= FORWARD_BACKWARD,
                                                             double saem_exponent= 1.,
                                                             int nb_iter= I_DEFAULT,
                                                             bool force_param= false) const;

/*   Hidden_markov_out_tree* hidden_markov_out_tree_estimation2(Format_error& error,
                                                              std::ostream& os,
                                                              const Hidden_markov_out_tree& ihmarkov,
                                                              bool counting_flag= true,
                                                              int state_trees= VITERBI,
                                                              int algorithm= FORWARD_BACKWARD,
                                                              double saem_exponent= 1.,
                                                              int nb_iter= I_DEFAULT,
                                                              bool force_param= false) const;*/

   /** Estimate a Hidden Markov Out Tree from the number of states and the data in \e self */
   Hidden_markov_out_tree* hidden_markov_out_tree_estimation(Format_error& error,
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

   void create_observation_histogram(int nb_state);
   // allocation of the histograms corresponding to the observation distributions
   void observation_histogram_computation();
   // computation of the histograms corresponding to the observation distributions
   void build_observation_histogram();
   // allocation and computation of the histograms
   // corresponding to the observation distributions

   // access to class members

   Hidden_markov_tree* get_markov() const; // { return markov; }
   Chain_data* get_chain_data() const; // { return chain_data; }
   double get_likelihood() const; // { return likelihood; }
   double get_hidden_likelihood() const; // { return hidden_likelihood; }
   int get_nb_states() const;

   ptHistogram_array_2d get_observation() const;
   ptHistogram_array get_observation(int variable) const;
   Histogram* get_observation(int variable, int state) const;

   Tree_characteristics* get_state_characteristics() const;
   // access to the characteristic quantity distributions
   // for the state variable

   ptOne_int_tree_set get_state_trees() const;
   Typed_edge_one_int_tree* get_state_tree(int itree) const;
   /** Return a Hidden_markov_tree_data containing the states
       as a variable */
   Hidden_markov_tree_data* get_state_hidden_markov_tree_data() const;
   /** For a given tree, and using a given type of entropy profile,
       return a Hidden_markov_tree_data containing the states,
       the entropy and the smoothed probabilities as variables */
   Hidden_markov_tree_data*
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
