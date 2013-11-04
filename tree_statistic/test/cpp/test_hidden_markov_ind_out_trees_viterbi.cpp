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
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
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
   typedef HiddenMarkovTreeData::state_tree_type state_tree_type;
   typedef HiddenMarkovTreeData::tree_type tree_type;
   typedef HiddenMarkovTreeData::key key;
   typedef HiddenMarkovTreeData::value value;
   typedef HiddenMarkovTree::double_array_2d double_array_2d;
   typedef HiddenMarkovTree::double_array_3d double_array_3d;
   typedef tree_type::vertex_descriptor vd;

   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   register unsigned int t, cptm= 0, nb_trees= 3;
   const int tid= 1, // id of the tree to be analyzed
             size= 55, // size2= 25;
             nb_children_max= 2;
   key vid;
   double hidden_likelihood, likelihood,
          likelihood2; //, max_marginal_entropy, entropy;
   value default_value;
   visitor v;
   vertex_array va;
   StatError error;
   const char * hmotpath= "./hmot_np_2s.hmt";// "./hmot.hmt";

   HiddenMarkovIndOutTree *hmot;
   HiddenMarkovTreeData *hmtd=NULL, *hmtdv= NULL, *smoothed= NULL,
                           *vud= NULL, *generalized= NULL,
                           *hmtdman= NULL, *hmtdtmp= NULL;
   state_tree_type **ptrees;
   tree_type **potrees;
   std::vector<ostringstream*> messages;
   std::vector<HiddenMarkovTreeData*> tree_list;

   // read and print a hidden Markov out tree
   hmot = hidden_markov_ind_out_tree_ascii_read(error, hmotpath);
   cout << error;

   if (hmot != NULL)
   {
      default_value.reset(0, 1);
      // default_value.Double(0)= .0;

      hmot->ascii_write(cout, false);
      cout << endl;

      // simulate a hidden Markov out tree
      hmtdtmp= hmot->simulation(error, nb_trees, size, nb_children_max);
      cout << error;
      // add a tree of size 2
      hmtdman= hmot->simulation(error, 1, 2, 1);
      cout << error;
      tree_list.push_back(hmtdman);
      hmtd= hmtdtmp->merge(error, tree_list);
      cout << error;
      nb_trees= hmtd->get_nb_trees();
      ptrees= new state_tree_type*[nb_trees];
      potrees= new tree_type*[nb_trees];

      // print the simulated trees
      cout << "Simulated trees : " << endl;
      for(t= 0; t < nb_trees; t++)
      {
         potrees[t]= hmtd->get_tree(t);
         potrees[t]->display(cout, 0);
         delete potrees[t];
      }
      cout << endl;
      delete [] potrees;
      potrees= NULL;

      cout << endl << "Simulated hidden trees : " << endl;
      for (t= 0; t < nb_trees; t++)
      {
         ptrees[t]= hmtd->get_state_tree(t);
         ptrees[t]->display(cout, 0);
         delete ptrees[t];
      }
      cout << endl;

      likelihood= hmot->likelihood_computation(*hmtd);
      hidden_likelihood= hmot->get_viterbi(*hmtd);

      // print the restored trees
      cout << endl << "Restored hidden trees (Viterbi): " << endl;
      for (t= 0; t < nb_trees; t++)
      {
         ptrees[t]= hmtd->get_state_tree(t);
         ptrees[t]->display(cout, 0);
         delete ptrees[t];
      }
      cout << endl;

      cout << "Completed likelihood: " << hidden_likelihood << endl;
      cout << "Likelihood: " << likelihood << endl;


      // Test of state_tree_computation
      hmtdv= hmot->state_tree_computation(error, *hmtd, VITERBI, true);

      cout << error << endl;
      if (hmtdv != NULL)
      {
         cout << endl << "Computation of the state trees: " << endl;
         for (t= 0; t < nb_trees; t++)
         {
            ptrees[t]= hmtdv->get_state_tree(t);
            ptrees[t]->display(cout, 0);
            delete ptrees[t];
         }

         cout << endl << "Restored hidden trees by smoothing: " << endl;
         // state restoration by smoothing
         hidden_likelihood= hmot->get_upward_downward(*hmtd);
         cout << endl;
         // print the restored trees
         for (t= 0; t < nb_trees; t++)
         {
            ptrees[t]= hmtd->get_state_tree(t);
            ptrees[t]->display(cout, 0);
            delete ptrees[t];
         }
         cout << endl;

         cout << "Completed likelihood : " << hidden_likelihood << endl;
         cout << "Likelihood : " << likelihood << endl;
         delete hmtdv;
         hmtdv= NULL;

#        if 0
         cout << endl << "Restored hidden trees by stochastic estimation: " << endl;
         hmtdv= hmot->sstate_simulation(*hmtd, hidden_likelihood);
         cout << endl;
         // print the restored trees
         if (hmtdv != NULL)
         {
            for (t= 0; t < nb_trees; t++)
            {
               ptrees[t]= hmtd->get_state_tree(t);
               ptrees[t]->display(cout, 0);
               delete ptrees[t];
            }
            cout << endl;

            cout << "Completed likelihood : " << hidden_likelihood << endl;
            delete hmtdv;
            hmtdv= NULL;
         }

         cout << endl << "Restored hidden trees by Gibbs sampling: " << endl;
         hmtdv= hmot->gibbs_state_simulation(*hmtd, hidden_likelihood);
         cout << endl;
         // print the restored trees
         if (hmtdv != NULL)
         {
            for (t= 0; t < nb_trees; t++)
            {
               ptrees[t]= hmtd->get_state_tree(t);
               ptrees[t]->display(cout, 0);
               delete ptrees[t];
            }
            cout << endl;

            cout << "Completed likelihood : " << hidden_likelihood << endl;
            delete hmtdv;
            hmtdv= NULL;
         }
#        endif
      }

      delete [] ptrees;
      ptrees= NULL;

      // Test of state_profile: UPWARD
      hmot->state_profile(error, *hmtd, tid, smoothed, hmtdv,
                          vud, generalized, messages,
                          GENERALIZED_VITERBI, 5,
                          UPWARD);
      cout << error << endl;
      ptrees= new state_tree_type*[1];

      if (smoothed != NULL)
      {
         cout << endl << "Computation of the 1st smoothed tree: (upward)"
              << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= smoothed->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete smoothed;
         smoothed= NULL;
      }
      if (vud != NULL)
      {
         cout << endl << "Viterbi upward-downward algorithm: " << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= vud->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete vud;
         vud= NULL;
         delete [] potrees;
         potrees= NULL;
      }

      if (generalized != NULL)
      {
         cout << endl << "Generalized viterbi: " << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= generalized->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete generalized;
         generalized= NULL;
         delete [] potrees;
         potrees= NULL;
      }

      // Test of state_profile: DOWNWARD
      hmot->state_profile(error, *hmtd, tid, smoothed, hmtdv,
                          vud, generalized, messages,
                          GENERALIZED_VITERBI, 2,
                          DOWNWARD);
      cout << error << endl;
      ptrees= new state_tree_type*[1];

      if (smoothed != NULL)
      {
         cout << endl << "Computation of the 1st smoothed tree (downward): "
              << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= smoothed->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete smoothed;
         smoothed= NULL;
      }
      if (vud != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete vud;
         vud= NULL;
      }

      if (generalized != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm]= NULL;
         delete generalized;
         generalized= NULL;
      }

      // Test of the generalized Viterbi algorithm
      // on subtrees
      vid= 2;
      generalized= hmot->generalized_viterbi_subtree(*hmtd, messages, 5,
                                                     likelihood, 0, vid);
      if (generalized != NULL)
      {
         cout << endl
              << "Generalized viterbi started at vertex "
              << vid << ": " << endl;
         cout << messages[++cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= generalized->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         /*cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;*/
         delete generalized;
         generalized= NULL;
         delete [] potrees;
         potrees= NULL;
      }

      // Test of the generalized Viterbi algorithm
      // on subtrees
      vid= 1;
      hmot->state_profile(error, *hmtd, tid, smoothed, hmtdv,
                          vud, generalized, messages,
                          GENERALIZED_VITERBI, 5,
                          UPWARD, vid);
      cout << error << endl;
      ptrees= new state_tree_type*[1];

      if (smoothed != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete smoothed;
         smoothed= NULL;
      }
      if (vud != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete vud;
         vud= NULL;
      }

      if (generalized != NULL)
      {
         cout << endl
              << "Generalized viterbi started at vertex "
              << vid << ": " << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= generalized->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete generalized;
         generalized= NULL;
         delete [] potrees;
         potrees= NULL;
      }
      // Test of the generalized Viterbi algorithm
      // on subtrees

      vid= 2;
      hmot->state_profile(error, *hmtd, tid, smoothed, hmtdv,
                          vud, generalized, messages,
                          GENERALIZED_VITERBI, 5,
                          UPWARD, vid);
      cout << error << endl;
      ptrees= new state_tree_type*[1];

      if (smoothed != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete smoothed;
         smoothed= NULL;
      }
      if (vud != NULL)
      {
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete vud;
         vud= NULL;
      }

      if (generalized != NULL)
      {
         cout << endl
              << "Generalized viterbi started at vertex "
              << vid << ": " << endl;
         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         potrees= new tree_type*[1];
         potrees[0]= generalized->get_tree(0);
         potrees[0]->display(cout, 0);
         cout << endl;
         delete potrees[0];
         potrees[0]= NULL;

         cout << messages[cptm]->str();
         delete messages[cptm];
         messages[cptm++]= NULL;
         delete generalized;
         generalized= NULL;
         delete [] potrees;
         potrees= NULL;
      }

      delete hmot;
      delete hmtd;
   }
   return 0;
}
