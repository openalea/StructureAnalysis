/****************************************************************
 *
 *  Test of the data structure and algorithms for hidden Markov
 *  out-trees as defined in hidden_markov_out_tree.h
 */

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/discrete_mixture.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/typed_edge_trees.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/int_trees.h"
#include "tree_statistic/hidden_markov_tree.h"
#include "tree_statistic/hidden_markov_ind_out_tree.h"

using namespace Stat_trees;

int main(void)
{
   typedef HiddenMarkovTreeData::tree_type tree_type;
   typedef HiddenMarkovTreeData::state_tree_type state_tree_type;
   typedef HiddenMarkovTreeData::value value;
   typedef HiddenMarkovTreeData::vertex_iterator vertex_iterator;
   typedef Stat_trees::HiddenMarkovTree::double_array_3d double_array_3d;
   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   const int nb_trees  = 5, size  = 25, nb_children_max  = 2;
   register int t; //, u, j, i, nb_states;
   value default_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   StatError error;
   const char *hmotpath  = "./hmot_np.hmt";
   const char *hmotparampath = "./hmot.hmt";
   int * const variables  = new int[1];
   HiddenMarkovIndOutTree *hmot= NULL, *hmot_estim= NULL;
   HiddenMarkovTree *ehmot = NULL;
   HiddenMarkovTreeData *hmtd, *state_hmtd, *ehmtd, *state_tree = NULL,
                        *hmtdv = NULL;
   tree_type **potrees  = NULL;
   state_tree_type **pstrees  = NULL;
   const DiscreteDistributionData *pdd  = NULL; //, *marginal=NULL;
   Trees *tmp_trees  = NULL;
   DiscreteMixtureData *mixture_data  = NULL;
   // double likelihood;

   // default constructor of HiddenMarkovTreeData
   hmtd = new HiddenMarkovTreeData();

   // destructor of HiddenMarkovTreeData
   delete hmtd;
   hmtd = NULL;

   potrees = new tree_type*[nb_trees];
   pstrees = new state_tree_type*[nb_trees];

   // read and print a non-parametric hidden Markov out tree
   cout << "Reference hidden Markov tree (non-parametric): " << endl;
   hmot= hidden_markov_ind_out_tree_ascii_read(error, hmotpath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   // simulation of hidden Markov out trees
   cout << "Simulate the attributes..." << endl;
   hmtd = hmot->simulation(error, nb_trees, size, nb_children_max);
   cout << error;

   cout << "Print the data with the hidden Markov tree: " << endl;
   cout << *hmtd << endl;

   hmot_estim= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot,
                                                           true, VITERBI, FORWARD_BACKWARD,
                                                           1., 10, true);
   cout << endl << "Parameter estimation:" << endl;
   cout << error;
   hmot_estim->ascii_write(cout, false);
   cout << "Print the data with the hidden Markov tree: " << endl;
   hmtd->ascii_write(cout, false);

   // print the simulated trees

   cout << endl << "Simulated trees : " << endl;
   for(t = 0; t < nb_trees; t++)
      if (t == 0)
      {
         potrees[t] = hmtd->get_tree(t);
         potrees[t]->display(cout, 0);
         delete potrees[t];
         potrees[t] = NULL;
      }

   cout << endl << "Simulated hidden trees : " << endl;
   for (t = 0; t < nb_trees; t++)
      if (t == 0)
      {
         pstrees[t] = hmtd->get_state_tree(t);
         pstrees[t]->display(cout, 0);
         delete pstrees[t];
         pstrees[t] = NULL;
      }
   cout << endl;

   cout << "Likelihood: " << hmtd->get_likelihood() << endl;

   cout << "Completed likelihood: " << hmtd->get_hidden_likelihood() << endl;

   cout << "Number of hidden states: " << hmtd->get_nb_states() << endl;

   cout << endl << "Marginal output histogram (1st variable) : " << endl;
   (hmtd->get_characteristics(0))->ascii_write_marginal_distribution(cout, true, false);

   cout << endl << "Marginal state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_marginal_distribution(cout, true, false);

   cout << endl << "Occurrence output histogram (1st variable): " << endl;
   (hmtd->get_characteristics(0))->ascii_write_nb_occurrences(cout, true, false);

   cout << endl << "Occurrence state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_nb_occurrences(cout, true, false);

   cout << endl << "State-conditional output histograms: " << endl;
   hmtd->ascii_write_observation(cout, true, false);

   cout << endl << "HiddenMarkovTreeData augmented with the state tree"
        <<" and the smoothed probabilities:" << endl;
   state_hmtd = hmtd->get_state_smoothed_hidden_markov_tree_data();
   state_hmtd->ascii_write(cout, false);

   cout << endl << "Extract marginal distribution of variable 1:" << endl;
   mixture_data = hmtd->extract_marginal(error, 1);
   cout << error << endl;
   if (mixture_data != NULL)
   {
      mixture_data->ascii_write(cout, false);
      delete mixture_data;
      mixture_data = NULL;
   }
   // Same thing using state_tree_computation
   hmtdv = hmtd->get_state_hidden_markov_tree_data();
   cout << endl
        << "Extract marginal distribution of variable 1 on restored tree:"
        << endl;
   mixture_data = hmtdv->extract_marginal(error, 2);
   cout << error << endl;
   if (mixture_data != NULL)
   {
      mixture_data->ascii_write(cout, false);
      delete mixture_data;
      mixture_data = NULL;
   }

   // extract markov part of hmtdv
   ehmot = hmtdv->extract_model(error);
   cout << error;
   if (ehmot != NULL)
   {
      cout << "Extract Markov part of HiddenMarkovTreeData:" << endl;
      ehmot->ascii_write(cout);
      delete ehmot;
      ehmot = NULL;
   }

   delete hmtdv;
   hmtdv = NULL;


#ifdef DEBUG
   cout << endl << "Plot of the augmented HiddenMarkovTreeData:" << endl;
   state_hmtd->plot_write(error, "/home/durand/tmpf");
   cout << error;
#endif

   ehmtd = hmot_estim->extract_data(error);
   cout << error;
   cout << endl << "Plain HiddenMarkovTreeData:" << endl;
   ehmtd->ascii_write(cout, false);
   cout << endl << "Extract an observation distribution from "
        << "plain HiddenMarkovTreeData:" << endl;
   pdd = ehmtd->extract(error, Stat_trees::OBSERVATION, 1, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract a observation distribution from "
        << "augmented HiddenMarkovTreeData:" << endl;
   pdd = state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "plain HiddenMarkovTreeData:" << endl;
   pdd = ehmtd->extract(error, Stat_trees::NB_ZONES, 1, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "augmented HiddenMarkovTreeData:" << endl;
   pdd = state_hmtd->extract(error, Stat_trees::NB_ZONES, 2, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   delete hmot;
   delete hmot_estim;
   delete hmtd;
   delete state_hmtd;

   hmot = NULL;
   hmot_estim = NULL;
   hmtd = NULL;
   state_hmtd = NULL;

   // read and print a parametric hidden Markov out tree
   cout << "Reference hidden Markov tree (parametric): " << endl;
   hmot= hidden_markov_ind_out_tree_ascii_read(error, hmotparampath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   // simulation of hidden Markov out trees
   cout << "Simulate the attributes..." << endl;
   hmtd = hmot->simulation(error, nb_trees, size, nb_children_max);
   cout << error;

   cout << "Print the data with the hidden Markov tree: " << endl;
   cout << *hmtd << endl;

   // extract markov part of state_hmtd
   ehmot = hmtd->extract_model(error);
   cout << error;
   if (ehmot != NULL)
   {
      cout << "Extract Markov part of HiddenMarkovTreeData:" << endl;
      ehmot->ascii_write(cout, true);
      delete ehmot;
      ehmot = NULL;
   }

   hmot_estim= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot,
                                                       true, VITERBI, FORWARD_BACKWARD,
                                                       1., 10, true);
   cout << endl << "Parameter estimation:" << endl;
   cout << error;
   hmot_estim->ascii_write(cout, false);
   cout << "Print the data with the hidden Markov tree: " << endl;
   hmtd->ascii_write(cout, false);

   // print the simulated trees

   cout << endl << "Simulated trees : " << endl;
   for(t = 0; t < nb_trees; t++)
      if (t == 0)
      {
         potrees[t] = hmtd->get_tree(t);
         potrees[t]->display(cout, 0);
         delete potrees[t];
         potrees[t] = NULL;
      }

   cout << endl << "Simulated hidden trees : " << endl;
   for (t = 0; t < nb_trees; t++)
      if (t == 0)
      {
         pstrees[t] = hmtd->get_state_tree(t);
         pstrees[t]->display(cout, 0);
         delete pstrees[t];
         pstrees[t] = NULL;
      }
   cout << endl;

   cout << "Likelihood: " << hmtd->get_likelihood() << endl;

   cout << "Completed likelihood: " << hmtd->get_hidden_likelihood() << endl;

   cout << "Number of hidden states: " << hmtd->get_nb_states() << endl;

   cout << endl << "Marginal output histogram (1st variable) : " << endl;
   (hmtd->get_characteristics(0))->ascii_write_marginal_distribution(cout, false, false);

   cout << endl << "Marginal state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_marginal_distribution(cout, false, false);

   cout << endl << "Occurrence output histogram (1st variable): " << endl;
   (hmtd->get_characteristics(0))->ascii_write_nb_occurrences(cout, false, false);

   cout << endl << "Occurrence state histogram: " << endl;
   (hmtd->get_state_characteristics())->ascii_write_nb_occurrences(cout, false, false);

   cout << endl << "State-conditional output histograms: " << endl;
   hmtd->ascii_write_observation(cout, false, false);

   cout << endl << "HiddenMarkovTreeData augmented with the state tree"
        <<" and the smoothed probabilities:" << endl;
   state_hmtd = hmtd->get_state_smoothed_hidden_markov_tree_data();
   state_hmtd->ascii_write(cout, false);
#ifdef DEBUG
   cout << endl << "Plot of the augmented HiddenMarkovTreeData:" << endl;
   state_hmtd->plot_write(error, "/home/durand/tmpg");
   cout << error;
#endif

   ehmtd = hmot_estim->extract_data(error);
   cout << error;
   cout << endl << "Plain HiddenMarkovTreeData:" << endl;
   ehmtd->ascii_write(cout, false);
   cout << endl << "Extract an observation distribution from "
        << "plain HiddenMarkovTreeData:" << endl;
   pdd = ehmtd->extract(error, Stat_trees::OBSERVATION, 1, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract an observation distribution from "
        << "augmented HiddenMarkovTreeData:" << endl;
   pdd = state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract an observation distribution from "
        << "augmented HiddenMarkovTreeData with a "
        << "wrong state identifier:" << endl;
   pdd = state_hmtd->extract(error, Stat_trees::OBSERVATION, 2, hmtd->get_nb_states());
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "plain HiddenMarkovTreeData:" << endl;
   pdd = ehmtd->extract(error, Stat_trees::NB_ZONES, 1, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   cout << endl << "Extract a characteristic distribution (zones) from "
        << "augmented HiddenMarkovTreeData:" << endl;
   pdd = state_hmtd->extract(error, Stat_trees::NB_ZONES, 2, 0);
   cout << error;
   if (pdd != NULL)
   {
      pdd->ascii_write(cout);
      delete pdd;
      pdd = NULL;
   }

   // estimation of an HMT from augmented HiddenMarkovTreeData
   delete hmot_estim;
   hmot_estim=NULL;

   variables[0]=1;
   tmp_trees = state_hmtd->select_variable(error, 1, variables);
   state_tree = new HiddenMarkovTreeData(*tmp_trees);
   delete tmp_trees;
   tmp_trees = NULL;
   hmot_estim= state_tree->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, true,
                                                             true, VITERBI, FORWARD_BACKWARD,
                                                             1., 0.99999, 1);
   delete [] variables;
   cout << error;
   hmot_estim->ascii_write(cout, false);

   delete state_tree;
   state_tree = NULL;
   delete hmot;
   delete hmot_estim;
   delete hmtd;
   delete state_hmtd;
   delete [] potrees;
   delete [] pstrees;

   hmot = NULL;
   hmot_estim = NULL;
   hmtd = NULL;
   state_hmtd = NULL;
   potrees = NULL;
   pstrees = NULL;
   return 0;
}
