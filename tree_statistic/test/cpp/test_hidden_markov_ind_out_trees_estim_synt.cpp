/****************************************************************
 *
 *  Test of the estimation algorithms for hidden Markov
 *  out-trees as defined in hidden_markov_ind_out_tree.h:
 *  checking the syntax and the algorithm basic properties
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
   typedef HiddenMarkovTreeData::tree_type hmtd_tree_type;
   typedef HiddenMarkovTreeData::pt_observed_trees_array pt_observed_trees_array;
   typedef HiddenMarkovTreeData::observed_trees observed_trees;
   typedef HiddenMarkovTreeData::value value;
   typedef HiddenMarkovTreeData::vertex_iterator vertex_iterator;
   typedef Stat_trees::HiddenMarkovTree::double_array_3d double_array_3d;
   typedef Stat_trees::HiddenMarkovTree::double_array_4d double_array_4d;
   typedef generic_visitor<hmtd_tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   int t;
   register int j, i, nb_states, nb_trees= 5;
   bool plot_ok= false;
   const int size= 15, nb_children_max= 2,
             nb_iterations= 5;
   unsigned int u;
   double likelihood, entropy1, entropy2, entropy3, entropyt,
          max_marginal_entropy;
   value default_value, current_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   StatError error;
   bool *force_param= NULL;
   int *perm;
   const char * hmotpath= "./hmot_np_2s.hmt";
   const char * hmotinitpath= "./hmot_np_init_2s.hmt";
   const char * hmotparampath= "./hmot_init.hmt";
   HiddenMarkovIndOutTree *hmot= NULL, *hmot2= NULL, *hmot_init= NULL,
                          *hmot_estim= NULL, *hmot_param= NULL;
   HiddenMarkovTreeData *hmtd= NULL, *hmtdman= NULL,
                           *hmtdtmp= NULL;
   Trees *ctrees= NULL;
   hmtd_tree_type **ptrees= NULL;
   observed_trees *merged_trees= NULL;
   std::vector<HiddenMarkovTreeData*> tree_list;
   hmtd_tree_type *ctree= NULL;
   Typed_edge_one_int_tree *citree= NULL;
   Unlabelled_typed_edge_tree *cstree= NULL;
   // const char * hsmcpath= "./laricio_3.hsc";
   double_array_3d state_marginal= NULL, output_cond= NULL,
          upward_prob= NULL, upward_parent_prob= NULL,
          downward_prob= NULL, state_entropy= NULL, state_entropy2= NULL,
          conditional_downward_entropy2 = NULL;
   double_array_4d downward_pair_prob= NULL, conditional_prob = NULL;
   double **sum_state_marginal= NULL, **conditional_upward_entropy = NULL,
          **conditional_downward_entropy = NULL, *upward_partial_entropy = NULL,
          **downward_partial_entropy1 = NULL, **downward_partial_entropy2 = NULL,
          *marginal_entropy = NULL;
   // int *rejected_variables= NULL;

   // reading and printing a hidden Markov out tree
   hmot= hidden_markov_ind_out_tree_ascii_read(error, hmotpath);
   hmot_param= hidden_markov_ind_out_tree_ascii_read(error, hmotparampath);

   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   if (hmot != NULL)
   {
      // simulation of hidden Markov out trees
      hmtdtmp= hmot->simulation(error, nb_trees, size, nb_children_max);
      cout << error;
      // add a tree of size 2
      hmtdman= hmot->simulation(error, 1, 2, 1);
      cout << error;
      tree_list.push_back(hmtdman);
      hmtd= hmtdtmp->merge(error, tree_list);
      delete hmtdtmp;
      hmtdtmp= NULL;
      delete hmtdman;
      hmtdman= NULL;

      // print indices for 1st tree
      cout << "Indices for tree number 0:" << endl;
      ctree= hmtd->get_tree(0);
      va= v.get_breadthorder(*ctree);
      for(u= 0; u < va.size(); u++)
      {
         current_value= ctree->get(va[u]);
         current_value.Int(0)= va[u];
         ctree->put(va[u], current_value);
      }
      ctree->display(cout, 0);
      delete ctree;
      ctree= NULL;

      if (hmtd != NULL)
      {
         nb_trees= hmtd->get_nb_trees();
         // print simulated trees

         cout << "Simulated trees : " << endl;
         for(t= 0; t < nb_trees; t++)
         {
            ctree= hmtd->get_tree(t);
            ctree->display(cout, 0);
            delete ctree;
            ctree= NULL;
         }

         cout << endl << "Simulated hidden trees : " << endl;
         for (t= 0; t < nb_trees; t++)
         {
            citree= hmtd->get_state_tree(t);
            citree->display(cout, 0);
            delete citree;
         }
         cout << endl;

         // compute the state marginal distributions
         hmot->get_state_marginal_distribution(*hmtd, state_marginal);
         sum_state_marginal= new double*[nb_trees];

         hmot->get_output_conditional_distribution(*hmtd, output_cond);

         nb_states= hmot->get_nb_state();


         ptrees= new hmtd_tree_type*[nb_trees];
         default_value.reset(0, 1);
         // default_value.Double(0)= .0;

         for(t= 0; t < nb_trees; t++)
         {
            if (t == 0)
               cout << "Marginal distributions for tree number " << t+1  << " : " << endl;
            citree= hmtd->get_state_tree(t);
            cstree= citree->get_structure();
            ptrees[t]= new hmtd_tree_type(*cstree, default_value);
            delete cstree;
            delete citree;
            va= v.get_breadthorder(*(ptrees[t]));
            sum_state_marginal[t]= new double[va.size()];
            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
               for(u= 0; u < va.size(); u++)
               {
                  default_value.Double(0)= state_marginal[t][j][va[u]];
                  ptrees[t]->put(va[u], default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
            if (t == 0)
               cout << "Output conditional distributions for tree number " << t+1  << " : " << endl;

            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
               for(u= 0; u < ptrees[t]->get_size(); u++)
               {
                  default_value.Double(0)= output_cond[t][j][u];
                  ptrees[t]->put(u, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());

            }
            for(u= 0; u < va.size(); u++)
            {
               sum_state_marginal[t][u]= .0;
               for(j= 0; j < hmot->get_nb_state(); j++)
                  sum_state_marginal[t][u]+= state_marginal[t][j][u];
            }
         }

         likelihood= hmot->get_upward_step(*hmtd, upward_prob, upward_parent_prob,
                                           state_entropy, state_marginal,
                                           output_cond, entropy1);
         hmot->get_downward_step(*hmtd, downward_prob, downward_pair_prob,
                                 upward_prob, upward_parent_prob,
                                 state_marginal, output_cond, entropy2);

         cout << "Log-likelihood of the parameter : " << likelihood << endl;
         // allocate conditional_downward_entropy2

         for(t= 0; t < nb_trees; t++)
         {
            entropyt= hmot->get_upward_partial_entropy(*hmtd, t, downward_prob,
                                                       state_entropy, upward_partial_entropy);

            if (t == 0)
               cout << "State conditional distributions for tree number " << t+1  << " : " << endl;
            for(j= 0; j < hmot->get_nb_state(); j++)
            {
               if (t == 0)
                  cout << "state " << j << " : " << endl;
           Tree_tie::tie(it, end)= ptrees[t]->vertices();
               while (it < end)
               {
                  default_value.Double(0)= downward_prob[t][j][*it];
                  ptrees[t]->put(*it++, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
            cout << "State tree entropy: " << entropyt << endl;
            cout << "Sum of conditional downward entropies: "
                 <<  hmot->get_downward_conditional_entropy(*hmtd, state_marginal, downward_prob,
                                                            upward_prob, upward_parent_prob,
                                                            conditional_downward_entropy,
                                                            conditional_downward_entropy2,
                                                            conditional_prob, state_entropy2, t) << endl;

            hmot->get_downward_partial_entropy(*hmtd, t, output_cond,
                                               state_marginal, upward_parent_prob,
                                               downward_prob, downward_pair_prob,
                                               state_entropy, conditional_downward_entropy2,
                                               conditional_prob,
                                               downward_partial_entropy2, false);

            cout << "Sum of conditional upward entropies: "
                 <<  hmot->get_upward_conditional_entropy(*hmtd, state_marginal, upward_prob,
                     upward_parent_prob, downward_prob,
                     conditional_upward_entropy, t) << endl;
            cout << "Sum of marginal entropies: "
                 << hmot->get_marginal_entropy(*hmtd, t, downward_prob,
                                               marginal_entropy, max_marginal_entropy) << endl;
            entropy3 = hmot->get_fast_downward_partial_entropy(*hmtd, t, downward_prob,
                                                               downward_pair_prob,
                                                               downward_partial_entropy1);

            cout << "State tree entropy (downward): " << entropy3 << endl;
            if (!hmot->downward_partial_entropy_comparison(*hmtd, t, output_cond,
                                                           state_marginal, upward_prob,
                                                           upward_parent_prob, downward_prob,
                                                           downward_pair_prob))
               cout << endl << "Raised a warning!" << endl;


            delete [] upward_partial_entropy;
            upward_partial_entropy = NULL;
            delete [] marginal_entropy;
            marginal_entropy = NULL;
         }

         default_value.reset(0, nb_states);
         for(t= 0; t < nb_trees; t++)
         // should be done for t < nb_trees
         {
            if (t == 0)
               cout << "Conditional distributions for pairs of states for tree number " << t+1  << " : " << endl;
            for(j= 0; j < nb_states; j++)
            {
               if (t == 0)
                  cout << "child state " << j << " : " << endl;
               Tree_tie::tie(it, end)= ptrees[t]->vertices();
               while (it < end)
               {
                  for(i= 0; i < nb_states; i++)
                     default_value.Double(i)= downward_pair_prob[t][i][j][*it];
                  ptrees[t]->put(*it++, default_value);
               }
               if (t == 0)
                  ptrees[t]->display(cout, ptrees[t]->root());
            }
         }

         delete hmot2;
         cout << "Testing estimation with initial hmt specification." << endl;
         cout << "Initialization HMT " << endl;
         hmot_init= hidden_markov_ind_out_tree_ascii_read(error, hmotinitpath);
         cout << error;
         if (hmot_init != NULL)
            hmot_init->ascii_write(cout, false);

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot_init,
                                                        true, VITERBI, FORWARD_BACKWARD,
                                                        1., 1, true);
         cout << error;

         if (hmot2 != NULL)
         {
            likelihood= hmot2->likelihood_computation(*hmtd);
            delete hmot_init;
            hmot_init= NULL;
            cout << "estimated HMT (one iteration)" << endl;
            // cout << *hmot2 << endl;
            hmot2->ascii_write(cout, false); //true
            for(i= 0; i < nb_iterations; i++)
            {
               hmot_init= hmot2;
               hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot_init,
                                                                  true, VITERBI,
                                                                  FORWARD_BACKWARD,
                                                                  1., 1, true);
               delete hmot_init;
               hmot_init= NULL;
               cout << "relative likelihood increase at iteration " << i << ": "
                    << (hmot2->likelihood_computation(*hmtd) - likelihood)/abs(likelihood) << endl;
               likelihood= hmot2->likelihood_computation(*hmtd);
            }
            cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
            cout << *hmot2;

            cout << endl << "Testing state profile plot" << endl;
            plot_ok= hmot2->state_profile_plot_write(error, "/home/durand/ftmp",
                                                     4, 8, NULL, UPWARD);
            plot_ok= plot_ok
               && hmot2->state_profile_plot_write(error, "/home/durand/ftmp",
                                                  4, 8, NULL, DOWNWARD);

            // permutation of the states
            perm= new int[hmot2->get_nb_state()];
            for(j= 0; j < hmot2->get_nb_state(); j++)
               perm[hmot2->get_nb_state()-j -1]=j;
            hmot2->state_permutation(error, perm);
            cout << error;
            cout << "Permutation of the states: " << endl;
            hmot2->ascii_write(cout, false);
            delete [] perm;
            perm= NULL;

            if (!plot_ok)
               cout << error;
            delete hmot2;
            hmot2= NULL;

         }
         cout << endl;

         cout << "Testing estimation with initial self-transition probability "
              << "specification." << endl;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, FORWARD_BACKWARD,
                                                        1., 0.99999, 1);
         cout << error;

         if (hmot2 != NULL)
         {
            likelihood= hmot2->likelihood_computation(*hmtd);
            hmot_init= NULL;
            cout << "estimated HMT (one iteration)" << endl;
            cout << *hmot2;
            for(i= 0; i < nb_iterations; i++)
            {
               hmot_init= hmot2;
               hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot_init,
                                                              true, VITERBI,
                                                              FORWARD_BACKWARD,
                                                              1., 1, true);
               delete hmot_init;
               hmot_init= NULL;
               cout << "relative likelihood increase at iteration " << i << ": "
                    << (hmot2->likelihood_computation(*hmtd) - likelihood)/abs(likelihood)<< endl;
               likelihood= hmot2->likelihood_computation(*hmtd);
            }
            cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
            cout << *hmot2;
            delete hmot2;
            hmot2= NULL;
         }

         cout << endl << "Testing estimation with auto-parametric distributions "
              << "and self-transitions." << endl;
         if (hmot_param != NULL)
         {
            // simulation of hidden Markov out trees
            hmtdtmp= hmot_param->simulation(error, nb_trees, size, nb_children_max);
            cout << error;
            if (hmtdtmp != NULL)
            {
                hmot2= hmtdtmp->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                                  true, VITERBI, FORWARD_BACKWARD,
                                                                  1., 0.99999, nb_iterations);
                cout << error;
                if (hmot2 != NULL)
                {
                   cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
                   cout << *hmot2;
                   delete hmot2;
                   hmot2= NULL;
                }
                delete hmtdtmp;
                hmtdtmp= NULL;
                delete hmot_param;
                hmot_param= NULL;
            }
         }


         cout << endl << "Testing estimation with imposed parametric distributions "
              << " and self-transitions." << endl;

         force_param= new bool[1];
         force_param[0]= true;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, FORWARD_BACKWARD,
                                                        1., 0.99999, nb_iterations, force_param);
         cout << error;
         delete [] force_param;
         force_param= NULL;

         if (hmot2 != NULL)
         {
            cout << endl << "estimated HMT (" << nb_iterations << " iterations)" << endl;
            cout << *hmot2;
            delete hmot2;
            hmot2= NULL;
         }

         cout << "Testing estimation using SEM algorithm" << endl;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, FORWARD_BACKWARD_SAMPLING,
                                                        0., 0.99999, 3);
         cout << error;
         if (hmot2 != NULL)
         {
            cout << *hmot2 << endl;
            delete hmot2;
            hmot2= NULL;
         }

         cout << "Testing estimation using EM algorithm 'a la Gibbs'" << endl;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, GIBBS_SAMPLING,
                                                        0., 0.99999, 3);
         cout << error;
         if (hmot2 != NULL)
         {
            cout << *hmot2 << endl;
            delete hmot2;
            hmot2= NULL;
         }

         cout << "Testing estimation using CEM algorithm" << endl;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                        true, VITERBI, VITERBI,
                                                        0., 0.99999, 3);
         cout << error;
         if (hmot2 != NULL)
         {
            cout << *hmot2 << endl;
            delete hmot2;
            hmot2= NULL;
         }

         cout << "Testing estimation using SAEM algorithm" << endl;

         hmot_estim= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', 2, false,
                                                             true, VITERBI, FORWARD_BACKWARD_SAMPLING,
                                                             0.5, 0.99999, 10);
         cout << error;
         if (hmot_estim != NULL)
         {
            cout << *hmot_estim << endl;
         }
         cout << endl;

         cout << "Testing estimation using SAEM algorithm "
              << "and initial model specification" << endl;

         hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot_estim,
                                                        true, VITERBI, FORWARD_BACKWARD_SAMPLING,
                                                        0.5, 10, true);
         cout << error;
         if (hmot2 != NULL)
         {
            cout << *hmot2 << endl;
            delete hmot2;
            hmot2= NULL;
         }
         if (hmot_estim != NULL)
         {
            delete hmot_estim;
            hmot_estim= NULL;
         }
         cout << endl;

         for(t= 0; t < nb_trees; t++)
         {
            for(j= 0; j < nb_states; j++)
            {
               delete [] output_cond[t][j];
               delete [] state_marginal[t][j];
               delete [] upward_prob[t][j];
               delete [] upward_parent_prob[t][j];
               delete [] downward_prob[t][j];
               for(i= 0; i < nb_states; i++)
               {
                  delete [] downward_pair_prob[t][j][i];
                  delete [] conditional_prob[t][j][i];
               }
               delete [] downward_pair_prob[t][j];
               delete [] conditional_prob[t][j];
               delete [] conditional_downward_entropy2[t][j];
               delete [] state_entropy[t][j];
               delete [] state_entropy2[t][j];
            }
            for(u= 0; u < ptrees[t]->get_size(); u++)
               if (fabs(sum_state_marginal[t][u] - 1.) > DOUBLE_ERROR)
                  cout << "State marginal distribution for tree " << t <<
                       " and node " << u << " sums to " << sum_state_marginal[t][u] << endl;
            delete [] state_marginal[t];
            delete [] output_cond[t];
            delete [] upward_prob[t];
            delete [] upward_parent_prob[t];
            delete [] downward_prob[t];
            delete [] downward_pair_prob[t];
            delete [] sum_state_marginal[t];
            delete [] conditional_upward_entropy[t];
            delete [] conditional_downward_entropy[t];
            delete [] conditional_downward_entropy2[t];
            delete [] conditional_prob[t];
            delete [] downward_partial_entropy1[t];
            delete [] downward_partial_entropy2[t];
            delete [] state_entropy[t];
            delete [] state_entropy2[t];
            delete ptrees[t];
         }

         delete [] ptrees;
         delete [] state_marginal;
         delete [] sum_state_marginal;
         delete [] output_cond;
         delete [] upward_prob;
         delete [] upward_parent_prob;
         delete [] downward_prob;
         delete [] downward_pair_prob;
         delete [] state_entropy;
         delete [] state_entropy2;
         delete [] conditional_prob;
         delete [] conditional_upward_entropy;
         delete [] conditional_downward_entropy;
         delete [] conditional_downward_entropy2;
         delete [] downward_partial_entropy1;
         delete [] downward_partial_entropy2;
         delete hmtd;
      }
      delete hmot;
   }
   return 0;
}
