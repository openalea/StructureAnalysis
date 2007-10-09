/****************************************************************
 *
 *  Comparison of a hidden Markov chain and a hidden Markov line-tree
 */

#include "tree_simple.h"
#include "tree_traits.h"
#include "basic_visitors.h"
#include "stat_tools.h"
#include "curves.h"
#include "markovian.h"
#include "sequences.h"
#include "regression.h" // used by markov.h
#include "markov.h"
#include "semi_markov.h"
#include "hidden_semi_markov.h"
#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "int_trees.h"
#include "hidden_markov_tree.h"
#include "hidden_markov_out_tree.h"

using namespace Stat_trees;

int main(void)
{
   typedef Hidden_markov_tree_data::state_tree_type state_tree_type;
   typedef Hidden_markov_tree_data::tree_type tree_type;
   typedef Hidden_markov_tree_data::value value;
   typedef Hidden_markov_tree_data::vertex_iterator vertex_iterator;
   typedef Stat_trees::Hidden_markov_tree::double_array_3d double_array_3d;
   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   register int var, t, u, nb_integral, cptm= 0; // j, i, nb_states,
   bool status;
   const int tid= 0; // id of the tree to be analyzed
   Hidden_semi_markov *hmc= NULL, *hmc_init= NULL;
   // Hidden_variable_order_markov *hmc= NULL, *hmc_init= NULL;
   Hidden_markov_out_tree *hmot= NULL, *hmot_init= NULL, *hmotref= NULL;
   Hidden_markov_tree_data *hmtd= NULL, *hmtdv= NULL, *smoothed= NULL,
                           *vud= NULL, *generalized= NULL;
   Markovian_sequences *ms= NULL;
   Sequences *seq= NULL;
   std::vector<ostringstream*> messages;
   tree_type *ctree, **potrees= NULL;
   state_tree_type **ptrees= NULL;
   Format_error error;
   const char * hmotrefpath= "./hmot_np.hmt";
   const char * hmotinitpath= "./hmot_np_init.hmt";
   const char * hmcinitpath= "./hmc_init.hmc";
   // const char * hmcinitpath= "./hmc_init.hvom";
   const int nb_trees= 10, size= 100, nb_children_max= 1;
   int ***sequences, *length;
   double likelihood, hidden_likelihood,
          check_likelihood, check_hidden_likelihood;
   visitor v;
   vertex_array va;


   // reading and printing of a hidden Markov out tree
   hmotref= hidden_markov_out_tree_ascii_read(error, hmotrefpath);
   cout << error;

   if (hmotref != NULL)
   {
      cout << "Reference hidden Markov model : " << endl;
      hmotref->ascii_write(cout, false);

      cout << endl;


      // simulation of hidden Markov out trees
      hmtd= hmotref->simulation(error, nb_trees, size, nb_children_max);
      cout << error;

#     ifdef DDEBUG
      // printing of the simulated trees
      cout << "Simulated trees : " << endl;
      for(t= 0; t < nb_trees; t++)
      {
         ctree= hmtd->get_tree(t);
         ctree->display(cout, 0);
      }

      cout << endl << "Simulated hidden trees : " << endl;
      for (t= 0; t < nb_trees; t++)
         (hmtd->get_state_tree(t))->display(cout, 0);
      cout << endl;
#     endif

      sequences= new int**[nb_trees];
      length= new int[nb_trees];
      nb_integral= hmtd->get_nb_int();
      for(t= 0; t < nb_trees; t++)
      {
         sequences[t]= new int*[nb_integral];
         ctree= hmtd->get_tree(t);
         length[t]= ctree->get_size();
         va= v.get_breadthorder(*ctree);
         for(var= 0; var < nb_integral; var++)
            sequences[t][var]= new int[length[t]];

         for(u= 0; u < length[t]; u++)
         {
            for(var= 0; var < nb_integral; var++)
               sequences[t][var][va[u]]= (ctree->get(va[u])).Int(var);
         }
      }

      ms= new Markovian_sequences(nb_integral, nb_trees, length, sequences);
#     ifdef DEBUG
      cout << "Markovian sequences: " << endl;
      ms->ascii_write(cout, false);

      cout << endl << "Plot of the Markovian Sequences into homedir:" << endl;
      ms->plot_write(error, "/home/jbdurand/tmp_ms");
      cout << error;

      cout << endl << "Plot of the Trees into homedir:" << endl;
      hmtd->Trees::plot_write(error, "/home/jbdurand/tmp_t");
      cout << error;
#     endif

      seq= hmtd->build_sequences(error, true);
      cout << endl << error;
#     ifdef DEBUG
      if (seq != NULL)
      {
         cout << "Sequences: " << endl;
         seq->ascii_write(cout, false);
         cout << endl << "Plot of the Sequences into homedir:" << endl;
         seq->plot_write(error, "/home/jbdurand/tmp_s");
         cout << error;
      }
#     endif

      cout << "Initialization HMT : " << endl;
      hmot_init= hidden_markov_out_tree_ascii_read(error, hmotinitpath);
      cout << error;
      hmot_init->ascii_write(cout, false);
      cout << endl;

      cout << "Initialization HMC : " << endl;
      // hmc_init= hidden_variable_order_markov_ascii_read(error, hmcinitpath);
      hmc_init= hidden_semi_markov_ascii_read(error, hmcinitpath);
      cout << error;
      if (hmc_init != NULL)
      {
         hmc_init->ascii_write(cout, false);
         cout << endl;

         hmot= hmtd->hidden_markov_out_tree_estimation(error, cout, *hmot_init,
                                                       true, VITERBI, 100, true);
         hmtd= hmot->get_markov_data();
         cout << error;

         hmc= ms->hidden_semi_markov_estimation(error, cout, *hmc_init, COMPLETE_LIKELIHOOD,
                                                true, VITERBI, 100);
         // hmc= ms->hidden_variable_order_markov_estimation(error, cout, *hmc_init, true,
         //                                                  true, VITERBI, 100);
         cout << error;

         // check the likelihood computation
         likelihood= hmtd->get_likelihood();
         check_likelihood= hmot->likelihood_computation(*hmtd);
         if (abs(likelihood-check_likelihood) > DOUBLE_ERROR)
            cout << "Warning: likelihood differs from Hidden_markov_tree_data"
                 << " to Hidden_markov_tree::likelihood_computation" << endl;

         // check the completed likelihood computation
         hidden_likelihood= hmtd->get_hidden_likelihood();
         check_hidden_likelihood= hmot->get_viterbi(*hmtd);
         if (abs(hidden_likelihood-check_hidden_likelihood) > DOUBLE_ERROR)
            cout << "Warning: completed likelihood differs from Hidden_markov_tree_data"
                 << " to Hidden_markov_tree::viterbi" << endl;
         cout << "Estimated HMT : " << endl;
         hmot->ascii_write(cout, false);
         cout << "log-likelihood of the state trees: " << hidden_likelihood << endl;
         cout << "log-likelihood of the observed trees: " << likelihood << endl;
         cout << endl;

         cout << "Estimated HMC : " << endl;
         hmc->ascii_write(cout, false);
         cout << endl;

         // check the state profile computation
         status= hmc->state_profile_ascii_write(error, cout, 1, GENERALIZED_VITERBI, 5);
         if (!status)
            cout << error << endl;

         status= hmot->state_profile(error, *hmtd, tid, smoothed, hmtdv,
                                     vud, generalized, messages, GENERALIZED_VITERBI, 5);
         if (!status)
            cout << error << endl;

         ptrees= new state_tree_type*[1];

         if (smoothed != NULL)
         {
            cout << endl << "Computation of the 1st smoothed tree: " << endl;
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

         for(t= 0; t < nb_trees; t++)
         {
            for(var= 0; var < nb_integral; var++)
               delete [] sequences[t][var];
            delete [] sequences[t];
         }
         delete [] sequences;
         delete [] length;

         if (hmotref != NULL)
            delete hmotref;
         if (hmtd != NULL)
            delete hmtd;
         if (seq != NULL)
            delete seq;
      }
   }
   return 0;
}
