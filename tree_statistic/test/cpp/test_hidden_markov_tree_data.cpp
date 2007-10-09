/****************************************************************
 *
 *  Test of the data structure and algorithms for hidden Markov
 *  out-trees as defined in hidden_markov_out_tree.h
 */

#include "tree_simple.h"
#include "tree_traits.h"
#include "basic_visitors.h"
#include "stat_tools.h"
#include "curves.h"
#include "markovian.h"
#include "sequences.h"
#include "semi_markov.h"
#include "distribution.h"
#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "int_trees.h"
#include "hidden_markov_tree.h"
#include "hidden_markov_out_tree.h"

using namespace Stat_trees;

int main(void)
{
   typedef Hidden_markov_tree_data::tree_type tree_type;
   typedef Hidden_markov_tree_data::state_tree_type state_tree_type;
   typedef Hidden_markov_tree_data::value value;
   typedef Hidden_markov_tree_data::vertex_iterator vertex_iterator;
   typedef Stat_trees::Hidden_markov_tree::double_array_3d double_array_3d;
   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   const int nb_trees= 5, size= 25, nb_children_max= 2;
   register int t; //, u, j, i, nb_states;
   value default_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   Format_error error;
   const char *hmotpath= "./hmot_np.hmt";
   const char *hmotparampath= "./hmot.hmt";
   int * const variables= new int[1];
   Hidden_markov_out_tree *hmot= NULL, *hmot_estim= NULL;
   Hidden_markov_tree_data *hmtd, *state_hmtd, *ehmtd, *state_tree=NULL;
   tree_type **potrees= NULL;
   state_tree_type **pstrees= NULL;
   const Distribution_data *pdd= NULL; //, *marginal=NULL;
   Trees *tmp_trees= NULL;
   // double likelihood;

   // default constructor of Hidden_markov_tree_data
   hmtd= new Hidden_markov_tree_data();

   // destructor of Hidden_markov_tree_data
   delete hmtd;
   hmtd= NULL;

   potrees= new tree_type*[nb_trees];
   pstrees= new state_tree_type*[nb_trees];

   // read and print a non-parametric hidden Markov out tree
   cout << "Reference hidden Markov tree (non-parametric): " << endl;
   hmot= hidden_markov_out_tree_ascii_read(error, hmotpath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   // simulation of hidden Markov out trees
   cout << "Simulate the attributes..." << endl;
   hmtd= hmot->simulation(error, nb_trees, size, nb_children_max);
   cout << error;

   cout << "Print the data with the hidden Markov tree: " << endl;
   cout << *hmtd << endl;

   hmot_estim= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot,
                                                       true, VITERBI, 10, true);
   cout << endl << "Parameter estimation:" << endl;
   cout << error;
   hmot_estim->ascii_write(cout, false);
   cout << "Print the data with the hidden Markov tree: " << endl;
   hmtd->ascii_write(cout, false);

   // print the simulated trees

   cout << endl << "Simulated trees : " << endl;
   for(t= 0; t < nb_trees; t++)
      if (t == 0)
      {
         potrees[t]= hmtd->get_tree(t);
         potrees[t]->display(cout, 0);
         delete potrees[t];
         potrees[t]= NULL;
      }

   cout << endl << "Simulated hidden trees : " << endl;
   for (t= 0; t < nb_trees; t++)
      if (t == 0)
      {
         pstrees[t]= hmtd->get_state_tree(t);
         pstrees[t]->display(cout, 0);
         delete pstrees[t];
         pstrees[t]= NULL;
      }
   cout << endl;

   cout << "Likelihood: " << hmtd->get_likelihood() << endl;

   cout << "Completed likelihood: " << hmtd->get_hidden_likelihood() << endl;

   cout << "Number of hidden states: " << hmtd->get_nb_states() << endl;

   cout << endl << "Marginal output histogram (1st variable) : " << endl;
   (hmtd->get_characteristics(0))->ascii_write_marginal(cout, true, false);

   cout << endl << "Marginal state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_marginal(cout, true, false);

   cout << endl << "Occurrence output histogram (1st variable): " << endl;
   (hmtd->get_characteristics(0))->ascii_write_nb_occurrences(cout, true, false);

   cout << endl << "Occurrence state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_nb_occurrences(cout, true, false);

   cout << endl << "State-conditional output histograms: " << endl;
   hmtd->ascii_write_observation(cout, true, false);

   cout << endl << "Hidden_markov_tree_data augmented with the state tree"
        <<" and the smoothed probabilities:" << endl;
   state_hmtd= hmtd->get_state_smoothed_hidden_markov_tree_data();
   state_hmtd->ascii_write(cout, false);

   /*
   cout << endl << "Marginal distribution of variable 1:" << endl;
   marginal= hmtd->extract_marginal(error, 1);
   marginal->ascii_write(cout, false);
   delete marginal;
   marginal= NULL; */

#ifdef DEBUG
   cout << endl << "Plot of the augmented Hidden_markov_tree_data:" << endl;
   state_hmtd->plot_write(error, "/home/jbdurand/tmpf");
   cout << error;
#endif

   ehmtd= hmot_estim->extract_data(error);
   cout << error;
   cout << endl << "Plain Hidden_markov_tree_data:" << endl;
   ehmtd->ascii_write(cout, false);
   cout << endl << "Extract an observation distribution from "
        << "plain Hidden_markov_tree_data:" << endl;
   pdd= ehmtd->extract(error, Stat_trees::OBSERVATION, 1, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract a observation distribution from "
        << "augmented Hidden_markov_tree_data:" << endl;
   pdd= state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "plain Hidden_markov_tree_data:" << endl;
   pdd= ehmtd->extract(error, Stat_trees::NB_ZONES, 1, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "augmented Hidden_markov_tree_data:" << endl;
   pdd= state_hmtd->extract(error, Stat_trees::NB_ZONES, 2, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   delete hmot;
   delete hmot_estim;
   delete hmtd;
   delete state_hmtd;

   hmot= NULL;
   hmot_estim= NULL;
   hmtd= NULL;
   state_hmtd= NULL;

   // read and print a parametric hidden Markov out tree
   cout << "Reference hidden Markov tree (parametric): " << endl;
   hmot= hidden_markov_out_tree_ascii_read(error, hmotparampath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   // simulation of hidden Markov out trees
   cout << "Simulate the attributes..." << endl;
   hmtd= hmot->simulation(error, nb_trees, size, nb_children_max);
   cout << error;

   cout << "Print the data with the hidden Markov tree: " << endl;
   cout << *hmtd << endl;

   hmot_estim= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot,
                                                       true, VITERBI, 10, true);
   cout << endl << "Parameter estimation:" << endl;
   cout << error;
   hmot_estim->ascii_write(cout, false);
   cout << "Print the data with the hidden Markov tree: " << endl;
   hmtd->ascii_write(cout, false);

   // print the simulated trees

   cout << endl << "Simulated trees : " << endl;
   for(t= 0; t < nb_trees; t++)
      if (t == 0)
      {
         potrees[t]= hmtd->get_tree(t);
         potrees[t]->display(cout, 0);
         delete potrees[t];
         potrees[t]= NULL;
      }

   cout << endl << "Simulated hidden trees : " << endl;
   for (t= 0; t < nb_trees; t++)
      if (t == 0)
      {
         pstrees[t]= hmtd->get_state_tree(t);
         pstrees[t]->display(cout, 0);
         delete pstrees[t];
         pstrees[t]= NULL;
      }
   cout << endl;

   cout << "Likelihood: " << hmtd->get_likelihood() << endl;

   cout << "Completed likelihood: " << hmtd->get_hidden_likelihood() << endl;

   cout << "Number of hidden states: " << hmtd->get_nb_states() << endl;

   cout << endl << "Marginal output histogram (1st variable) : " << endl;
   (hmtd->get_characteristics(0))->ascii_write_marginal(cout, false, false);

   cout << endl << "Marginal state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_marginal(cout, false, false);

   cout << endl << "Occurrence output histogram (1st variable): " << endl;
   (hmtd->get_characteristics(0))->ascii_write_nb_occurrences(cout, false, false);

   cout << endl << "Occurrence state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_nb_occurrences(cout, false, false);

   cout << endl << "State-conditional output histograms: " << endl;
   hmtd->ascii_write_observation(cout, false, false);

   cout << endl << "Hidden_markov_tree_data augmented with the state tree"
        <<" and the smoothed probabilities:" << endl;
   state_hmtd= hmtd->get_state_smoothed_hidden_markov_tree_data();
   state_hmtd->ascii_write(cout, false);
#ifdef DEBUG
   cout << endl << "Plot of the augmented Hidden_markov_tree_data:" << endl;
   state_hmtd->plot_write(error, "/home/jbdurand/tmpg");
   cout << error;
#endif

   ehmtd= hmot_estim->extract_data(error);
   cout << error;
   cout << endl << "Plain Hidden_markov_tree_data:" << endl;
   ehmtd->ascii_write(cout, false);
   cout << endl << "Extract an observation distribution from "
        << "plain Hidden_markov_tree_data:" << endl;
   pdd= ehmtd->extract(error, Stat_trees::OBSERVATION, 1, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract an observation distribution from "
        << "augmented Hidden_markov_tree_data:" << endl;
   pdd= state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract an observation distribution from "
        << "augmented Hidden_markov_tree_data with a "
        << "wrong state identifier:" << endl;
   pdd= state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, hmtd->get_nb_states());
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "plain Hidden_markov_tree_data:" << endl;
   pdd= ehmtd->extract(error, Stat_trees::NB_ZONES, 1, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "augmented Hidden_markov_tree_data:" << endl;
   pdd= state_hmtd->extract(error, Stat_trees::NB_ZONES, 2, 0);
   cout << error;
   if (pdd!= NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd= NULL;
   }

   // estimation of an HMT from augmented Hidden_markov_tree_data
   delete hmot_estim;
   hmot_estim=NULL;

   variables[0]=1;
   tmp_trees= state_hmtd->select_variable(error, 1, variables);
   state_tree= new Hidden_markov_tree_data(*tmp_trees);
   delete tmp_trees;
   tmp_trees= NULL;
   hmot_estim= state_tree->hidden_markov_out_tree_estimation(error, cout, 'o', 2, true,
                                                             true, VITERBI, 0.99999, 1);
   delete [] variables;
   cout << error;
   hmot_estim->ascii_write(cout, false);

   delete state_tree;
   state_tree= NULL;
   delete hmot;
   delete hmot_estim;
   delete hmtd;
   delete state_hmtd;
   delete [] potrees;
   delete [] pstrees;

   hmot= NULL;
   hmot_estim= NULL;
   hmtd= NULL;
   state_hmtd= NULL;
   potrees= NULL;
   pstrees= NULL;
   return 0;
}
