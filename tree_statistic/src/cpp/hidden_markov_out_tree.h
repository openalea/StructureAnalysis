/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/STAT_TREES/src/hidden_markov_out_tree.h,v $
 *       $Id: hidden_markov_out_tree.h 2722 2007-02-16 14:17:56Z jbdurand $
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

#ifndef HIDDEN_MARKOV_OUT_TREE_H
#define HIDDEN_MARKOV_OUT_TREE_H

namespace Stat_trees
{

/*! \file hidden_markov_out_tree.h
    \brief Purpose:
     provide the class "Hidden_markov_out_tree", which corresponds
     to hidden Markov trees with independent children states,
     given their parent
*/

/****************************************************************
 *
 *  Class definitions :
 */

/**
   \class Hidden_markov_out_tree
   \brief a hidden Markov out-tree model, with independent
    children states, given their parent
*/
class Hidden_markov_out_tree : public Hidden_markov_tree
{  // a hidden Markov out-tree model, with independent
   // children states given their parent

   friend class Hidden_markov_tree_data;

   friend Hidden_markov_out_tree* hidden_markov_out_tree_ascii_read(Format_error& error,
                                                                    const char * path,
                                                                    int size,
                                                                    bool counting_flag,
                                                                    double cumul_threshold);

   friend ostream& operator<<(ostream &os, const Hidden_markov_out_tree& markov);

protected :

   ostream& ascii_write(ostream& os, const Hidden_markov_tree_data* otrees,
                        bool exhaustive= false, bool file_flag= false,
                        const Test* test= NULL) const;
   ostream& spreadsheet_write(ostream& os, const Hidden_markov_tree_data * tree,
                              const Test * test= NULL) const;
   bool plot_write(const char * prefix, const char * title,
                   const Hidden_markov_tree_data * otrees) const;

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

   void state_marginal_distribution(const Hidden_markov_tree_data& trees,
                                    double_array_3d& state_marginal,
                                    int index= I_DEFAULT) const;
   double** state_marginal_distribution(const Trees& trees,
                                        int index) const;
   /*
   void output_conditional_distribution(const Hidden_markov_tree_data& trees,
                                        double_array_3d& output_cond,
                                        bool logcomputation= false,
                                        int index= I_DEFAULT) const;*/
   /**Return the parameter likelihood,
      and compute the state tree process entropy*/
   double upward_step(const Hidden_markov_tree_data& trees,
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
   void downward_step(const Hidden_markov_tree_data& trees,
                      double_array_3d& downward_prob, double_array_4d& downward_pair_prob,
                      double_array_3d upward_prob, double_array_3d upward_parent_prob,
                      double_array_3d marginal_prob,
                      double_array_3d output_cond_prob,
                      double& entropy2,
                      int index= I_DEFAULT) const;
                      // should be const double *** const marginal_prob, etc

   double likelihood_correction(const Hidden_markov_tree_data& trees) const;

   /**Compute the maximal marginal entropy, the state tree entropy
   and return the completed likelihood */
   double upward_downward(const Hidden_markov_tree_data& trees,
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
   double smoothed_probabilities(const Hidden_markov_tree_data& trees,
                                 double_array_3d& smoothed_prob,
                                 double_array_2d& marginal_entropy,
                                 double_array_2d& conditional_entropy,
                                 double_array_2d& partial_entropy,
                                 int index,
                                 int entropy_algo= UPWARD) const;

   double viterbi(const Hidden_markov_tree_data& trees,
                  int index= I_DEFAULT) const;

   /**Compute the number of possible state trees */
   long double nb_state_trees(const Hidden_markov_tree_data& trees,
                              int index= I_DEFAULT) const;

   /**Upward_downward Viterbi algorithm */
   Hidden_markov_tree_data* viterbi_upward_downward(const Hidden_markov_tree_data& trees,
                                                    std::vector<ostringstream*>& messages,
                                                    double likelihood,
                                                    double& state_likelihood,
                                                    std::deque<int>*& vd,
                                                    int index= I_DEFAULT,
                                                    std::ostream* os= NULL,
                                                    char format= 'a',
                                                    int vertex= I_DEFAULT) const;

   /**Compute the nb_state_trees best state trees starting under some vertex*/

   Hidden_markov_tree_data* generalized_viterbi_subtree(const Hidden_markov_tree_data& trees,
                                                        std::vector<ostringstream*>& messages,
                                                        int nb_state_trees,
                                                        double likelihood,
                                                        int index,
                                                        int root) const;

   /**Compute the nb_state_trees best state trees */
   Hidden_markov_tree_data* generalized_viterbi(const Hidden_markov_tree_data& trees,
                                                std::vector<ostringstream*>& messages,
                                                int nb_state_trees,
                                                double likelihood,
                                                int index= I_DEFAULT,
                                                int sroot= I_DEFAULT) const;

   /**Simulate a state tree given observed tree (SEM-like principle)*/
   Hidden_markov_tree_data* sstate_simulation(const Hidden_markov_tree_data& trees,
                                              double& state_likelihood,
                                              bool characteristic_flag= true) const;

   /**Simulate a state tree given observed tree, via Gibbs sampling*/
   Hidden_markov_tree_data* gibbs_state_simulation(const Hidden_markov_tree_data& trees,
                                                   double& state_likelihood,
                                                   bool characteristic_flag= true) const;

   /**Compute the completed likelihood for a given set of trees and state trees*/
   double state_likelihood_computation(const Hidden_markov_tree_data& trees) const;

   /**Compute the completed likelihood for one tree and one state tree*/
   double state_likelihood_computation(const Hidden_markov_tree_data& trees,
                                       int index) const;

   /** Compute the marginal entropy */
   void marginal_entropy_computation(const Hidden_markov_tree_data& trees,
                                     int t,
                                     double_array_3d downward_prob,
                                     double*& marginal_entropy,
                                     double& max_marginal_entropy) const;

   /** Compute the entropy of partial state processes: partial state subtree
       rooted at each possible node*/
   double upward_partial_entropy_computation(const Hidden_markov_tree_data& trees,
                                             int t,
                                             double_array_3d downward_prob,
                                             double_array_3d state_entropy,
                                             double*& partial_entropy) const;

   /** Compute the entropy of partial state processes: state subtree
       pruned at subtrees rooted at each possible node*/
   double downward_partial_entropy_computation(const Hidden_markov_tree_data& trees,
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
   void upward_conditional_entropy_computation(const Hidden_markov_tree_data& trees,
                                               double_array_3d marginal_prob,
                                               double_array_3d upward_prob,
                                               double_array_3d upward_parent_prob,
                                               double_array_3d downward_prob,
                                               double_array_2d& conditional_entropy,
                                               int index= I_DEFAULT) const;

   /** Compute the conditional entropy: entropy of each state given
       all non-descendant states*/
   void downward_conditional_entropy_computation(const Hidden_markov_tree_data& trees,
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

   Hidden_markov_out_tree();
   /** Create a Hidden Markov Out Tree from the type of transition probability matrix,
   the number of states, the number of integer and floating observation processes,
   the number of values for each integer process
   and an indicator array for forcing parametric integer processes */
   Hidden_markov_out_tree(char itype, int inb_state, int inb_ioutput_process,
                          int inb_doutput_process, int* nb_value,
                          bool* force_param= NULL);
   Hidden_markov_out_tree(const Chain * pchain, int inb_ioutput_process,
                          Nonparametric_process** pobservation,
                          int size, bool counting_flag);
   Hidden_markov_out_tree(const Chain * pchain,
                          int inb_ioutput_process, int inb_doutput_process,
                          Nonparametric_process** nonparametric_observation,
                          Parametric_process** iparametric_observation,
                          Parametric_process** dparametric_observation,
                          int size, bool counting_flag);
   Hidden_markov_out_tree(const Hidden_markov_out_tree& markov, bool data_flag= true,
                          bool characteristic_flag= true);
   Hidden_markov_out_tree(const Hidden_markov_tree& markov, bool data_flag= true,
                          bool characteristic_flag= true);
   Hidden_markov_out_tree(const Hidden_markov_out_tree& markov, double self_transition);
   Hidden_markov_out_tree(const Hidden_markov_out_tree& markov, char manip, int param);
   virtual ~Hidden_markov_out_tree();

   /** Return the data part of a Hidden_markov_tree,
       keeping a reference on \e self */
   Hidden_markov_tree_data* extract_data(Format_error& error) const;
   ostream& line_write(ostream& os) const;
   ostream& ascii_write(ostream& os, bool exhaustive= false) const;
   bool ascii_write(Format_error& error, const char * path,
                    bool exhaustive= false) const;
   bool spreadsheet_write(Format_error& error, const char * path) const;
   bool plot_write(Format_error& error, const char * prefix,
                   const char * title= NULL) const;

   double likelihood_computation(const Trees& trees, int index) const;
   double likelihood_computation(const Hidden_markov_tree_data& trees) const;

   /** Compute optimal state trees */
   Hidden_markov_tree_data* state_tree_computation(Format_error& error,
                                                   const Trees& trees,
                                                   int algorithm= VITERBI,
                                                   bool characteristic_flag= false,
                                                   int index= I_DEFAULT,
                                                   int entropy_algo= UPWARD) const;

   /**Compute the state profiles for a given tree,
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
                      int entropy_algo= UPWARD,
                      int root= I_DEFAULT) const;

   /** Write Gnuplot output of state, entropy and Viterbi profiles */
   bool tree_state_profile_plot_write(Format_error &error, const char *prefix,
                                      const Hidden_markov_tree_data& trees,
                                      int identifier, int vertex,
                                      const char *title= NULL,
                                      int entropy_algo= UPWARD) const;

   Hidden_markov_tree_data* simulation(Format_error& error, const Histogram& ihsize,
                                       const Histogram& ihnb_children,
                                       bool counting_flag= true,
                                       // bool divergence_flag= true) const;
                                       bool divergence_flag= false) const;
   Hidden_markov_tree_data* simulation(Format_error& error, int inb_trees,
                                       int isize, int inb_children_max,
                                       bool counting_flag= true) const;
   Hidden_markov_tree_data* simulation(Format_error& error, int inb_trees,
                                       const Trees& otrees,
                                       bool counting_flag= true) const;

   // methods below should be used for testing purposes only...
   void get_state_marginal_distribution(const Hidden_markov_tree_data& trees,
                                        double_array_3d& state_marginal) const;
   void get_output_conditional_distribution(const Hidden_markov_tree_data& trees,
                                            double_array_3d& output_cond) const;
   double get_upward_step(const Hidden_markov_tree_data& trees, double_array_3d& upward_prob,
                          double_array_3d& upward_parent_prob, double_array_3d& state_entropy,
                          double_array_3d marginal_prob, double_array_3d output_cond_prob,
                          double& entropy1) const;
                          // should be const double *** const marginal_prob, etc
   void get_downward_step(const Hidden_markov_tree_data& trees,
                          double_array_3d& downward_prob, double_array_4d& downward_pair_prob,
                          double_array_3d upward_prob, double_array_3d upward_parent_prob,
                          double_array_3d marginal_prob,
                          double_array_3d output_cond_prob, double& entropy2) const;
   double get_viterbi(const Hidden_markov_tree_data& trees,
                      int index= I_DEFAULT);
   /* Return the completed likelihood */
   double get_upward_downward(const Hidden_markov_tree_data& trees, int index= I_DEFAULT,
                              ostream* os= NULL, char format= 'a') const;

};
   Hidden_markov_out_tree* hidden_markov_out_tree_ascii_read(Format_error& error,
                                                             const char * path,
                                                             int size= I_DEFAULT_TREE_SIZE,
                                                             bool counting_flag= true,
                                                             double cumul_threshold= OCCUPANCY_THRESHOLD);
};

#endif
