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
   typedef HiddenMarkovTreeData::tree_type tree_type;
   typedef HiddenMarkovTreeData::value value;
   typedef generic_visitor<tree_type> visitor;
   typedef visitor::vertex_array vertex_array;

   register int t, j;
   unsigned int u;
   const int nb_trees= 5, size= 15, nb_children_max= 2;
   bool w= false;
   value default_value;
   visitor v;
   vertex_array va;
   StatError error;
   const char * hsmcpath= "./laricio_3.hsc";
   // const char * hmtpath= "./hmt.hmt";
   const char * hmotpath= "./hmot.hmt";
   const char * hmotpath2= "./hidden_markov.hmt";
   const char * wapath= "./hmot_ascii.hmt";
   const char * wspath= "./hmot_spreadsheet.hmt";
   const char * swapath= "./laricio_3ascii.hsc";
   // const char * swspath= "./laricio_3spreadsheet.hsc";
   int* perm;
   double*** state_marginal= NULL, ***output_cond= NULL;
   double** sum_state_marginal= NULL;
   HiddenSemiMarkov *hsmc;
   // HiddenMarkovTree *hmt, *hmt2;
   HiddenMarkovIndOutTree *hmot, *hmot2;
   HiddenMarkovTreeData *hmtd;
   tree_type **ptrees;

   // default constructor of HiddenMarkovTreeData
   hmtd= new HiddenMarkovTreeData();

   // destructor of HiddenMarkovTreeData
   delete hmtd;
   hmtd= NULL;

   // read and print a hidden semi-Markov chain
   hsmc= hidden_semi_markov_ascii_read(error, hsmcpath);
   cout << error;

   if (hsmc != NULL)
   {
      hsmc->ascii_write(cout, false);
      cout << endl;

      w= hsmc->ascii_write(error, swapath);
      if (!w)
         cout << error;

      // does not work for the moment
      /*
      w= hsmc->spreadsheet_write(error, swspath);
      if (!w)
         cout << error; */

      // read and print a hidden Markov out tree (error ?)

      hmot = Stat_trees::hidden_markov_ind_out_tree_ascii_read(error, hmotpath2);
      cout << error;

      if (hmot != NULL)
      {
         // hmot->ascii_write(cout, false);
         delete hmot;
         hmot = NULL;
      }

      // read and print a hidden Markov out tree
      hmot= Stat_trees::hidden_markov_ind_out_tree_ascii_read(error, hmotpath);
      cout << error;
      cerr << error;

      if (hmot != NULL)
         hmot->ascii_write(cout, false);

      // permutation of the states
      perm= new int[hmot->get_nb_state()];
      for(j= 0; j < hmot->get_nb_state(); j++)
         perm[hmot->get_nb_state()-j -1]=j;
      hmot->state_permutation(error, perm);
      cout << error;
      cout << "Permutation of the states: " << endl;
      hmot->ascii_write(cout, false);
      delete [] perm;
      perm= NULL;

      cout << endl;

      if (hmot != NULL)
      {
         // Copy constructor of a hidden Markov out tree
         // without the characteristics
         hmot2= new HiddenMarkovIndOutTree(*hmot, false, false);

         hmot2->ascii_write(cout, false);
         cout << endl;

         delete hmot2;
         hmot2= NULL;
         // write a hidden Markov tree into a file :
         w= hmot->ascii_write(error, wapath);
         if (!w)
            cout << error;

         w= hmot->spreadsheet_write(error, wspath);
         if (!w)
            cout << error;

         // simulate a hidden Markov out tree
         hmtd= hmot->simulation(error, nb_trees, size, nb_children_max);
         cout << error;

         // print the simulated trees

         cout << "Simulated trees : " << endl;
         for(t= 0; t < nb_trees; t++)
            (hmtd->get_tree_ptr(t))->display(cout, 0);

         cout << endl << "Simulated hidden trees : " << endl;
         for (t= 0; t < nb_trees; t++)
            (hmtd->get_state_tree_ptr(t))->display(cout, 0);
         cout << endl;

         // compute the state marginal distributions
         hmot->get_state_marginal_distribution(*hmtd, state_marginal);
         sum_state_marginal= new double*[nb_trees];

         hmot->get_output_conditional_distribution(*hmtd, output_cond);

         ptrees= new tree_type*[nb_trees];
         default_value.reset(0, 1);
         // default_value.Double(0)= .0;

         for(t= 0; t < nb_trees; t++)
         // should be done for t < nb_trees
         {
            if (t == 0)
               cout << "Marginal distributions for tree number " << t+1  << " : " << endl;
            ptrees[t]= new tree_type(*((hmtd->get_state_tree(t))->get_structure()),
                                     default_value);
            va= v.get_breadthorder(*(ptrees[t]));
            /* cout << "Tree number " << t+1 << " is supposed to have " << va.size()
                 << " nodes." << endl; */
            sum_state_marginal[t]= new double[va.size()];
            sum_state_marginal[t][0]= .0;
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

         for(t= 0; t < nb_trees; t++)
         {
            for(j= 0; j < hmot->get_nb_state(); j++)
               delete [] state_marginal[t][j];
            for(u= 0; u < ptrees[t]->get_size(); u++)
               if (fabs(sum_state_marginal[t][u] - 1) > DOUBLE_ERROR)
                  cout << "State marginal distribution for tree " << t <<
                       " and node " << u << " sums to " << sum_state_marginal[t][u] << endl;
            delete [] state_marginal[t];
            delete [] sum_state_marginal[t];
            delete ptrees[t];
         }

         // Copy constructor of a hidden Markov out tree
         // with the characteristics
         hmot2= new HiddenMarkovIndOutTree(*hmot, true, true);

         /* delete hmt;
         delete hmt2; */
         delete [] ptrees;
         delete [] state_marginal;
         delete [] sum_state_marginal;
         delete hsmc;
         delete hmot2;
      }

      delete hmot;
      // delete hmtd;
      return 0;
   }
}
