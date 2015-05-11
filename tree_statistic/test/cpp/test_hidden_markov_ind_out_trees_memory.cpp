/****************************************************************
 *
 *  Test memory leaks in hidden Markov
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

   const int nb_trees = 1, size = pow(2,10)+2,  nb_children_max = 2,
             nb_simulation = 2, //40*40*100,
             v_size = pow(2, 26);
   register int t; //, u, j, i, nb_states;
   value default_value;
   visitor v;
   vertex_iterator it, end;
   vertex_array va;
   StatError error;
   const char *hmotpath  = "./hmot_np.hmt";
   const char *hmotparampath = "./hmot.hmt";
   HiddenMarkovIndOutTree *hmot= NULL, *hmot_estim= NULL;
   HiddenMarkovTree *ehmot = NULL;
   HiddenMarkovTreeData *hmtd, *state_hmtd, *ehmtd, *state_tree = NULL,
                        *hmtdv = NULL;
   tree_type **potrees  = NULL;
   state_tree_type **pstrees  = NULL;
   const DiscreteDistributionData *pdd  = NULL; //, *marginal=NULL;
   Trees *tmp_trees  = NULL;
   DiscreteMixtureData *mixture_data  = NULL;
   Int_fl_container ifc;

   // read and print an Int_fl_container
   /* cout << "Initialize Int_fl_container: " << endl;
   for(t = 0; t < 10; t++)
   {
      // if ((t % 10) == 0)
      cout << t << "\t";
      ifc = Int_fl_container(v_size, v_size);
   }
   cout << endl; */

   // read and print a non-parametric hidden Markov out tree
   cout << "Simulate hidden Markov tree (non-parametric): " << endl;
   hmot = hidden_markov_ind_out_tree_ascii_read(error, hmotpath);
   cout << error;

   if (hmot != NULL)
      hmot->ascii_write(cout, false);

   cout << endl;

   cout << "Simulate the attributes once" << endl;
   hmtd = hmot->simulation(error, nb_trees, size, nb_children_max);
   cout << error;
   if (hmtd != NULL)
   {
      hmtd->ascii_write(cout, false);
      delete hmtd;
      hmtd = NULL;
   }

   cout << "Simulate the attributes " << nb_simulation << " times: ";
   // simulation of hidden Markov out trees
   for(t = 0; t < nb_simulation; t++)
   {
      if ((t % 50) == 0)
         cout << t << "\t";
      hmtd = hmot->simulation(error, nb_trees, size, nb_children_max);
      delete hmtd;
      hmtd = NULL;
   }
   cout << endl << error;

   delete hmot;
   delete hmtd;

   hmot = NULL;
   hmtd = NULL;

   return 0;
}
