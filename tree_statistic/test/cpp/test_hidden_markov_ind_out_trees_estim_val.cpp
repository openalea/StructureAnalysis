/****************************************************************
 *
 *  Test of the estimation algorithms for hidden Markov
 *  out-trees as defined in hidden_markov_out_tree.h:
 *  checking the returned values
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
   typedef HiddenMarkovTreeData::tree_type tree_type;
   typedef HiddenMarkovTreeData::value value;
   typedef HiddenMarkovTreeData::vertex_iterator vertex_iterator;
   typedef Stat_trees::HiddenMarkovTree::double_array_3d double_array_3d;
   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   HiddenMarkovIndOutTree *hmot= NULL, *hmot2= NULL, *hmot_init= NULL;
   HiddenMarkovTreeData *hmtd;
   // tree_type **ptrees;
   tree_type ctree;
   value default_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   StatError error;
   // const char * hsmcpath= "./laricio_3.hsc";
   // const char * hmtpath= "./hmt.hmt";
   const char * hmotpath= "./hmot.hmt";
   const char * hmotinitpath= "./hmot_init.hmt";
   const int nb_children_max= 2, nb_trees= 50, size= 40,
             nb_iterations= 10;
   register int nb_states= 0; //t, u, j, i,
   // double_array_3d upward_prob= NULL, upward_parent_prob= NULL,
   //        downward_prob= NULL, *downward_pair_prob= NULL;
   //        state_marginal= NULL, output_cond= NULL,
   // double **sum_state_marginal= NULL;
   // double likelihood;

   // read and print a hidden Markov out tree
   hmot= hidden_markov_ind_out_tree_ascii_read(error, hmotpath);

   // Copy constructor of a hidden Markov out tree
   if (hmot != NULL)
   {
      cout << "Read a hidden Markov tree from file " << hmotpath << endl;
      hmot2= new HiddenMarkovIndOutTree(*hmot, false, false);

      hmot2->ascii_write(cout, false);
      cout << endl;

      cout << "Simulate trees from this hidden Markov tree" << endl;
      hmtd= hmot->simulation(error, size, nb_trees, nb_children_max);
      cout << error;

      cout << "Initialization HMT: " << endl;
      hmot_init= hidden_markov_ind_out_tree_ascii_read(error, hmotinitpath);
      cout << error;
      if (hmot_init != NULL)
         hmot_init->ascii_write(cout, false);

      hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, *hmot_init,
                                                     true, VITERBI,
                                                     FORWARD_BACKWARD,
                                                     1., nb_iterations, true);
      cout << error;

      if (hmot2 != NULL)
      {
         cout << endl << "Estimated HMT using above initialization:" << endl;
         hmot2->ascii_write(cout, false);
         nb_states= hmot2->get_nb_state();
         delete hmot2;
         hmot2=NULL;
      }

      if (hmot2 != NULL)
      {
         cout << endl << "Estimated HMT using above initialization (basic EM algorithm):" << endl;
         hmot2->ascii_write(cout, false);
         nb_states= hmot2->get_nb_state();
         delete hmot2;
         hmot2=NULL;
      }

      hmot2= hmtd->hidden_markov_ind_out_tree_estimation(error, cout, 'o', nb_states,
                                                     false, true, VITERBI,
                                                     FORWARD_BACKWARD,
                                                     1., 0.7, nb_iterations);
      cout << error;
      cout << endl;

      if (hmot2 != NULL)
      {
         cout << endl << "Estimated HMT using initial self-transitions:" << endl;
         hmot2->ascii_write(cout, false);
         delete hmot2;
         hmot2=NULL;
      }

      delete hmot;
      delete hmtd;
   }
   return 0;
}
