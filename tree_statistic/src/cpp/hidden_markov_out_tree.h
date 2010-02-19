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
 *       $Id: hidden_markov_out_tree.h 2722 2007-02-16 14:17:56Z jbdurand $
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

#ifndef HIDDEN_MARKOV_OUT_TREE_H
#define HIDDEN_MARKOV_OUT_TREE_H

namespace Stat_trees
{

/*! \file hidden_markov_out_tree.h
    \brief Purpose:
     provide the class "HiddenMarkovOutTree", which corresponds
     to hidden Markov trees with independent children states,
     given their parent
*/

/****************************************************************
 *
 *  Class definitions :
 */

/**
   \class HiddenMarkovOutTree
   \brief a hidden Markov out-tree model, with independent
    children states, given their parent
*/
class HiddenMarkovOutTree : public HiddenMarkovTree
{  // a hidden Markov out-tree model, with independent
   // children states given their parent

   friend class HiddenMarkovTreeData;

   friend HiddenMarkovOutTree* hidden_markov_out_tree_ascii_read(StatError& error,
                                                                 const char * path,
                                                                 int size,
                                                                 bool counting_flag,
                                                                 double cumul_threshold);

   friend ostream& operator<<(ostream &os, const HiddenMarkovOutTree& markov);

protected :

   ostream& ascii_write(ostream& os, const HiddenMarkovTreeData* otrees,
                        bool exhaustive= false, bool file_flag= false,
                        const Test* test= NULL) const;
   ostream& spreadsheet_write(ostream& os, const HiddenMarkovTreeData * tree,
                              const Test * test= NULL) const;
   bool plot_write(const char * prefix, const char * title,
                   const HiddenMarkovTreeData * otrees) const;

   void state_no_occurrence_probability(int state, double increment= LEAVE_INCREMENT);
   void state_first_occurrence_root_distribution(int state, int min_nb_value= 1,
                                                 double cumul_threshold= CUMUL_THRESHOLD);
   void state_first_occurrence_leaves_distribution(int state, int min_nb_value= 1,
                                                   double cumul_threshold= CUMUL_THRESHOLD);
   void state_leave_probability(const double * memory, int state,
                                double increment= LEAVE_INCREMENT);
   void state_sojourn_size_distribution(const double * memory, int state,
                                        int min_nb_value= 1,
                                        double cumul_threshold= CUMUL_THRESHOLD);
   void state_nb_pattern_mixture(int state, char pattern);

   void output_no_occurrence_probability(int variable, int output,
                                         double increment= LEAVE_INCREMENT);
   void output_first_occurrence_root_distribution(int variable, int output,
                                                  int min_nb_value= 1,
                                                  double cumul_threshold= CUMUL_THRESHOLD);
   void output_first_occurrence_leaves_distribution(int variable, int output,
                                                    int min_nb_value= 1,
                                                    double cumul_threshold= CUMUL_THRESHOLD);
   void output_leave_probability(const double * memory,
                                 int variable, int output,
                                 double increment= LEAVE_INCREMENT);
   void output_sojourn_size_distribution(const double * memory, int variable,
                                         int output, int min_nb_value= 1 ,
                                         double cumul_threshold= CUMUL_THRESHOLD);
   void output_nb_zones_mixture(int variable, int output);
   void output_nb_occurrences_mixture(int variable, int output);

   int nb_parameter_computation(double min_probability= 0.) const;

   void state_marginal_distribution(const HiddenMarkovTreeData& trees,
                                    double_array_3d& state_marginal,
                                    int index= I_DEFAULT) const;
   double** state_marginal_distribution(const Trees& trees,
                                        int index) const;
   /*
   void output_conditional_distribution(const HiddenMarkovTreeData& trees,
                                        double_array_3d& output_cond,
                                        bool logcomputation= false,
                                        int index= I_DEFAULT) const;*/
   /**Return the parameter likelihood,
      and compute the state tree process entropy*/
   double upward_step(const HiddenMarkovTreeData& trees,
                      double_array_3d& upward_prob,
                      double_array_3d& upward_parent_prob,
                      double_array_3d& state_entropy,
                      double_array_3d marginal_prob,
                      double_array_3d output_cond_prob,
                      double& entropy1,
                      int index= I_DEFAULT) const;
                      // should be const double *** const marginal_prob, etc

   /**Compute the state posterior probabilities,
      the state posterior probabilities,
      and the state tree process entropy - loglikelihood*/
   void downward_step(const HiddenMarkovTreeData& trees,
                      double_array_3d& downward_prob, double_array_4d& downward_pair_prob,
                      double_array_3d upward_prob, double_array_3d upward_parent_prob,
                      double_array_3d marginal_prob,
                      double_array_3d output_cond_prob,
                      double& entropy2,
                      int index= I_DEFAULT) const;
                      // should be const double *** const marginal_prob, etc

   double likelihood_correction(const HiddenMarkovTreeData& trees) const;

   /**Compute the maximal marginal entropy, the state tree entropy
   and return the completed likelihood */
   double upward_downward(const HiddenMarkovTreeData& trees,
                          double& max_marginal_entropy,
                          double& entropy1,
                          double& likelihood,
                          std::deque<int>*& vd,
                          int index= I_DEFAULT,
                          std::ostream* os= NULL,
                          char format= 'a',
                          int vertex= I_DEFAULT,
                          int entropy_algo= UPWARD) const;

   /**Compute the smoothed probabilities and return the likelihood */
   double smoothed_probabilities(const HiddenMarkovTreeData& trees,
                                 double_array_3d& smoothed_prob,
                                 double_array_2d& marginal_entropy,
                                 double_array_2d& conditional_entropy,
                                 double_array_2d& partial_entropy,
                                 int index,
                                 int entropy_algo= UPWARD) const;

   double viterbi(const HiddenMarkovTreeData& trees,
                  int index= I_DEFAULT) const;

   /**Compute the number of possible state trees */
   long double nb_state_trees(const HiddenMarkovTreeData& trees,
                              int index= I_DEFAULT) const;

   /**Upward_downward Viterbi algorithm */
   HiddenMarkovTreeData* viterbi_upward_downward(const HiddenMarkovTreeData& trees,
                                                 std::vector<ostringstream*>& messages,
                                                 double likelihood,
                                                 double& state_likelihood,
                                                 std::deque<int>*& vd,
                                                 int index= I_DEFAULT,
                                                 std::ostream* os= NULL,
                                                 char format= 'a',
                                                 int vertex= I_DEFAULT) const;

   /**Compute the nb_state_trees best state trees */
   HiddenMarkovTreeData* generalized_viterbi(const HiddenMarkovTreeData& trees,
                                             std::vector<ostringstream*>& messages,
                                             int nb_state_trees,
                                             double likelihood,
                                             int index= I_DEFAULT,
                                             int sroot= I_DEFAULT) const;

   /**Simulate a state tree given observed tree (SEM-like principle)*/
   HiddenMarkovTreeData* sstate_simulation(const HiddenMarkovTreeData& trees,
                                           double& state_likelihood,
                                           bool characteristic_flag= true) const;

   /**Simulate a state tree given observed tree, via Gibbs sampling*/
   HiddenMarkovTreeData* gibbs_state_simulation(const HiddenMarkovTreeData& trees,
                                                double& state_likelihood,
                                                bool characteristic_flag= true) const;

   /**Compute the completed likelihood for a given set of trees and state trees*/
   double state_likelihood_computation(const HiddenMarkovTreeData& trees) const;

   /**Compute the completed likelihood for one tree and one state tree*/
   double state_likelihood_computation(const HiddenMarkovTreeData& trees,
                                       int index) const;

   /** Compute the marginal entropy */
   void marginal_entropy_computation(const HiddenMarkovTreeData& trees,
                                     int t,
                                     double_array_3d downward_prob,
                                     double*& marginal_entropy,
                                     double& max_marginal_entropy) const;

   /** Compute the entropy of partial state processes: partial state subtree
       rooted at each possible node*/
   double upward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                             int t,
                                             double_array_3d downward_prob,
                                             double_array_3d state_entropy,
                                             double*& partial_entropy) const;

   /** Compute the entropy of partial state processes: state subtree
       pruned at subtrees rooted at each possible node*/
   double downward_partial_entropy_computation(const HiddenMarkovTreeData& trees,
                                               int t,
                                               double_array_3d output_cond_prob,
                                               double_array_3d marginal_prob,
                                               double_array_3d upward_parent_prob,
                                               double_array_3d downward_prob,
                                               double_array_3d state_entropy,
                                               double_array_3d conditional_entropy,
                                               double_array_4d conditional_prob,
                                               double*& partial_entropy) const;

   /** Compute the conditional entropy: entropy of each state given
       the children state subtrees*/
   void upward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                               double_array_3d marginal_prob,
                                               double_array_3d upward_prob,
                                               double_array_3d upward_parent_prob,
                                               double_array_3d downward_prob,
                                               double_array_2d& conditional_entropy,
                                               int index= I_DEFAULT) const;

   /** Compute the conditional entropy: entropy of each state given
       all non-descendant states*/
   void downward_conditional_entropy_computation(const HiddenMarkovTreeData& trees,
                                                 double_array_3d marginal_prob,
                                                 double_array_3d downward_prob,
                                                 double_array_3d upward_prob,
                                                 double_array_3d upward_parent_prob,
                                                 double_array_2d& expected_conditional_entropy,
                                                 double_array_3d& conditional_entropy,
                                                 double_array_4d& conditional_prob,
                                                 double_array_3d& state_entropy,
                                                 int index= I_DEFAULT) const;

public :

   HiddenMarkovOutTree();
   /** Create a Hidden Markov Out Tree from the type of transition probability matrix,
   the number of states, the number of integer and floating observation processes,
   the number of values for each integer process
   and an indicator array for forcing parametric integer processes */
   HiddenMarkovOutTree(char itype, int inb_state, int inb_ioutput_process,
                       int inb_doutput_process, int* nb_value,
                       bool* force_param= NULL);
   HiddenMarkovOutTree(const Chain * pchain, int inb_ioutput_process,
                       NonparametricProcess** pobservation,
                       int size, bool counting_flag);
   HiddenMarkovOutTree(const Chain * pchain,
                       int inb_ioutput_process, int inb_doutput_process,
                       NonparametricProcess** nonparametric_observation,
                       DiscreteParametricProcess** iparametric_observation,
                       DiscreteParametricProcess** dparametric_observation,
                       int size, bool counting_flag);
   HiddenMarkovOutTree(const HiddenMarkovOutTree& markov, bool data_flag= true,
                       bool characteristic_flag= true);
   HiddenMarkovOutTree(const HiddenMarkovTree& markov, bool data_flag= true,
                       bool characteristic_flag= true);
   HiddenMarkovOutTree(const HiddenMarkovOutTree& markov, double self_transition);
   HiddenMarkovOutTree(const HiddenMarkovOutTree& markov, char manip, int param);
   virtual ~HiddenMarkovOutTree();

   /** Return the data part of a HiddenMarkovTree,
       keeping a reference on \e self */
   HiddenMarkovTreeData* extract_data(StatError& error) const;
   ostream& line_write(ostream& os) const;
   ostream& ascii_write(ostream& os, bool exhaustive= false) const;
   bool ascii_write(StatError& error, const char * path,
                    bool exhaustive= false) const;
   bool spreadsheet_write(StatError& error, const char * path) const;
   bool plot_write(StatError& error, const char * prefix,
                   const char * title= NULL) const;

   double likelihood_computation(const Trees& trees, int index) const;
   double likelihood_computation(const HiddenMarkovTreeData& trees) const;

   /** Compute optimal state trees */
   HiddenMarkovTreeData* state_tree_computation(StatError& error,
                                                const Trees& trees,
                                                int algorithm= VITERBI,
                                                bool characteristic_flag= false,
                                                int index= I_DEFAULT,
                                                int entropy_algo= UPWARD) const;

   /** Compute the nb_state_trees best state trees starting under some vertex */

   HiddenMarkovTreeData* generalized_viterbi_subtree(const HiddenMarkovTreeData& trees,
                                                     std::vector<ostringstream*>& messages,
                                                     int nb_state_trees,
                                                     double likelihood,
                                                     int index,
                                                     int root) const;

   /**Compute the state profiles for a given tree,
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

   /** Write Gnuplot output of state, entropy and Viterbi profiles */
   bool tree_state_profile_plot_write(StatError &error, const char *prefix,
                                      const HiddenMarkovTreeData& trees,
                                      int identifier, int vertex,
                                      const char *title= NULL,
                                      int entropy_algo= UPWARD) const;

   HiddenMarkovTreeData* simulation(StatError& error, const FrequencyDistribution& ihsize,
                                    const FrequencyDistribution& ihnb_children,
                                    bool counting_flag= true,
                                    // bool divergence_flag= true) const;
                                    bool divergence_flag= false) const;
   HiddenMarkovTreeData* simulation(StatError& error, int inb_trees,
                                    int isize, int inb_children_max,
                                    bool counting_flag= true) const;
   HiddenMarkovTreeData* simulation(StatError& error, int inb_trees,
                                    const Trees& otrees,
                                    bool counting_flag= true) const;

   // methods below should be used for testing purposes only...
   void get_state_marginal_distribution(const HiddenMarkovTreeData& trees,
                                        double_array_3d& state_marginal) const;
   void get_output_conditional_distribution(const HiddenMarkovTreeData& trees,
                                            double_array_3d& output_cond) const;
   double get_upward_step(const HiddenMarkovTreeData& trees, double_array_3d& upward_prob,
                          double_array_3d& upward_parent_prob, double_array_3d& state_entropy,
                          double_array_3d marginal_prob, double_array_3d output_cond_prob,
                          double& entropy1) const;
                          // should be const double *** const marginal_prob, etc
   void get_downward_step(const HiddenMarkovTreeData& trees,
                          double_array_3d& downward_prob, double_array_4d& downward_pair_prob,
                          double_array_3d upward_prob, double_array_3d upward_parent_prob,
                          double_array_3d marginal_prob,
                          double_array_3d output_cond_prob, double& entropy2) const;
   double get_viterbi(const HiddenMarkovTreeData& trees,
                      int index= I_DEFAULT);
   /* Return the completed likelihood */
   double get_upward_downward(const HiddenMarkovTreeData& trees, int index= I_DEFAULT,
                              ostream* os= NULL, char format= 'a') const;

};
   HiddenMarkovOutTree* hidden_markov_out_tree_ascii_read(StatError& error,
                                                          const char * path,
                                                          int size= I_DEFAULT_TREE_SIZE,
                                                          bool counting_flag= true,
                                                          double cumul_threshold= OCCUPANCY_THRESHOLD);
};

#endif
