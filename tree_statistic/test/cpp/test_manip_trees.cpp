/****************************************************************
 *
 *  Test of the data structure for observed trees as defined
 *  in observed_trees.h
 */

#include <math.h>

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "sequence_analysis/sequences.h"
#include "tree_statistic/int_fl_containers.h"
#include "tree_statistic/typed_edge_trees.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include "tree_statistic/int_trees.h"

#ifdef _WIN32
#define random rand
#endif

using namespace Stat_trees;
using namespace Tree_tie;

int main(void)
{

   typedef Int_trees::tree_type tree_type;
   typedef Trees::tree_type rtree_type;
   typedef Int_trees::int_trees int_trees;
   typedef Int_trees::int_array int_array;
   typedef Int_trees::pt_int_trees_array pt_int_trees_array;
   typedef tree_type::vertex_descriptor vertex;
   typedef tree_type::pt_Distribution_array pt_Distribution_array;
   typedef tree_type::value value;
   typedef tree_type::vertex_iterator vertex_iterator;
   const int ainf_bound= 1, asup_bound= 20,
             sinf_bound= 1, ssup_bound= 4;
   const double probability= 0.6;
   const int ident= UNIFORM;
   const double parameter= D_DEFAULT;
   const int nb_trees= 1, nb_var_trees= 2;
   const DiscreteParametric sdistrib(ident, sinf_bound, ssup_bound, parameter, probability);
   const DiscreteParametric adistrib(ident, ainf_bound, asup_bound, parameter, probability);
   const int max_depth= 3;
   const int max_size= 20;
   register int t, o, *itype= NULL, *citype= NULL, i;
   const register int nb_class= 2;
   value v, rv;
   Unlabelled_typed_edge_tree *tmp_utree= NULL;
   tree_type *default_base_tree= NULL, tmp_base_tree;
   tree_type  **otrees= NULL;
   vertex_iterator it, end;
   Int_trees *b= NULL;
   rtree_type **ctrees= NULL;
   tree_type **citrees= NULL;
   Sequences *res_seq= NULL;
   Trees *c= NULL, *d= NULL;
   pt_int_trees_array pota;
   pt_Distribution_array pda;
   int_trees *res= NULL, *merged= NULL, *cmerged= NULL, *clustered= NULL;
   int_array limit, symbol_array, iidentifier, ivariables;
   Typed_edge_one_int_tree *ident_tree= NULL; //*otrees1= NULL,
   Trees** potc= new Trees*[0];
   tree_type ti;
   vertex *va= new vertex[10];
   StatError error;

   pota= new int_trees*[nb_var_trees];
   pda= new Distribution*[2];
   // We check that the basic structure work

   itype= new int[2];
   itype[0]= INT_VALUE;
   itype[1]= INT_VALUE;

   citype= new int[2];
   citype[0]= INT_VALUE;
   citype[1]= REAL_VALUE;

   // Random assignement of the structures

   pda[0]= new Distribution(adistrib);
   pda[1]= pda[0];

   v.reset(2,0);
   cout << endl << "Tree attributes : " << v << endl;
   // We check that the basic structures work

   cout << "Default tree : " << endl;
   default_base_tree= new tree_type(2, 0, 0, 1);
   default_base_tree->add_vertex(v);
   // necessary : default_base_tree was empty at previous stage
   default_base_tree->display(cout,default_base_tree->root());

   tmp_utree= new Unlabelled_typed_edge_tree;
   otrees= new tree_type*[nb_trees];

   for(t= 0; t < nb_trees; t++)
   {
      tmp_utree->simulation(sdistrib, max_size, max_depth);
      otrees[t]= new tree_type(*default_base_tree);
      otrees[t]->set_structure(*tmp_utree, v);
   }

   // build and simulate a reference forest
   b= new Int_trees(nb_trees, itype, otrees);

   b->iid_simulation(adistrib);

   // merge reference forest with a collection of forests
   for(o= 0; o < nb_var_trees; o++)
   {
      for(t= 0; t < nb_trees; t++)
      {
         tmp_utree->simulation(sdistrib, max_size, max_depth);
         delete otrees[t];
         otrees[t]= new tree_type(*default_base_tree);
         otrees[t]->set_structure(*tmp_utree, v);
         otrees[t]->iid_simulation(pda);
      }
      pota[o]= new Int_trees(nb_trees, itype, otrees);
   }

   for(t= 0; t < nb_trees; t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }
   delete [] otrees;
   otrees= NULL;

   otrees= new tree_type*[b->get_nb_trees()];

   cout << "Trees to be merged : " << endl;
   for(t= 0; t < b->get_nb_trees(); t++)
   {
      otrees[t]= b->get_tree(t);
      otrees[t]->display(cout, otrees[t]->root());
      cout << endl;
   }
   for(i= 0; i < 2; i++)
      for(t= 0; t < pota[i]->get_nb_trees(); t++)
      {
         delete otrees[t];
         otrees[t]= NULL;
         otrees[t]= pota[i]->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
      }

   cout << "Merge trees..." << endl;
   merged= b->merge(error, nb_var_trees, pota);
   cmerged= new int_trees(*merged);
   cout << error;

   for(t= 0; t < b->get_nb_trees(); t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }

   delete [] otrees;
   otrees= NULL;

   if (merged != NULL)
   {
      cout << endl << "Merged trees: " << endl;
      otrees= new tree_type*[merged->get_nb_trees()];
      merged->ascii_write(cout, false);
      for(t= 0; t < merged->get_nb_trees(); t++)
      {
         otrees[t]= merged->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
      }
      delete [] otrees;
      otrees= NULL;
   }

   // shift reference forest
   cout << "Shift trees..." << endl;

   res= b->shift(error, 1, -1);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Shifted tree: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         delete otrees[t];
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   // cluster reference forest using a step
   cout << "Cluster trees..." << endl;

   res= b->cluster(error, 2, 5);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Tree with 2nd variable clustered: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         delete otrees[t];
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   limit= new int[nb_class];
   limit[0]= 6;
   limit[nb_class-1]= asup_bound;

   // cluster reference forest using an array of limits
   cout << "Cluster trees..." << endl;
   res= b->cluster(error, 1, nb_class, limit);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Tree with 1st variable clustered "
           << "( " << nb_class << " clusters with bounds = [0, "
           << 6 << ", " <<  asup_bound << "]): " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete [] otrees;
      otrees= NULL;
   }

   symbol_array= new int[nb_class];
   symbol_array[0]= 3;
   symbol_array[1]= 2;

   // transcode reference forest
   cout << "Transcode trees..." << endl;
   res= res->transcode(error, 1, symbol_array);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Tree with 1st variable transcoded (0->3, "
           << "1->2): " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   // select some trees within merged forest
   // based on the values
   cout << "Select value..." << endl;
   res= merged->value_select(error, 1, 0, 2, true);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Trees with 1st variable between 0 and 2: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   // select some trees within merged forest
   // based on the identifiers
   iidentifier= new int[nb_class];
   iidentifier[0]= 0;
   iidentifier[1]= 2;
   cout << "Select individual..." << endl;
   res= merged->select_individual(error, 2, iidentifier, false);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Trees with identifier different from 0 or 2: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   // select variables in merged forest
   ivariables= new int[1];
   ivariables[0]= 2;
   // ivariables[1]= 1;
   cout << "Select variable..." << endl;
   res= merged->select_variable(error, 1, ivariables, true);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Selection of the 2nd variable" << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         delete otrees[t];
         otrees[t]= NULL;
         cout << endl;
      }
      delete [] otrees;
      otrees= NULL;
      // cluster values within selected variables
      cout << endl << "Clustering with step 2 on the selected variable" << endl;
      clustered= res->cluster(error, 1, 2);
      if (clustered != NULL)
      {
         otrees= new tree_type*[clustered->get_nb_trees()];
         clustered->ascii_write(cout, false);
         for(t= 0; t < res->get_nb_trees(); t++)
         {
            otrees[t]= clustered->get_tree(t);
            otrees[t]->display(cout, otrees[t]->root());
            delete otrees[t];
            otrees[t]= NULL;
            cout << endl;
         }
         delete [] otrees;
         otrees= NULL;

         for(o= 0; o < nb_var_trees; o++)
         {
            delete pota[o];
            pota[o]= NULL;
         }
         delete [] pota;

         pota= new int_trees*[1];
         pota[0]= res;

         // merge variables of the forest with clustered variables and
         // that of the forest with selected variables
         cout << "Merge variable..." << endl;
         res= clustered->merge_variable(error, 1, pota);
         cout << error;

         if (res != NULL)
         {
            // cout << endl << "Merging twice the variable" << endl;
            otrees= new tree_type*[res->get_nb_trees()];
            res->ascii_write(cout, false);
            for(t= 0; t < res->get_nb_trees(); t++)
            {
               otrees[t]= res->get_tree(t);
               otrees[t]->display(cout, otrees[t]->root());
               delete otrees[t];
               otrees[t]= NULL;
            }
            delete res;
            res= NULL;
            delete [] otrees;
            otrees= NULL;
         }
         delete clustered;
         clustered= NULL;
      }
   }

   cout << "Select size..." << endl;
   res= merged->size_select(error, 10, 25);
   cout << error;

   if (res != NULL)
   {
      cout << endl << "Trees with size between 10 and 25: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   cout << "Get tree identifiers..." << endl;
   ident_tree= merged->get_identifier_tree(0);

   if (ident_tree != NULL)
   {
      cout << endl << "Tree identifiers for 1st tree: " << endl;
      ident_tree->display(cout, ident_tree->root());
      delete ident_tree;
      ident_tree= NULL;
   }

   cout << "Select subtrees..." << endl;
   res= merged->select_subtrees(error, 0, 2, false);
   cout << error;

   delete merged;
   merged= NULL;

   if (res != NULL)
   {
      cout << endl << "Selection of subtrees rooted at node 3: " << endl;
      otrees= new tree_type*[res->get_nb_trees()];
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete res;
      res= NULL;
      delete [] otrees;
      otrees= NULL;
   }

   // test of the StatError messages
   // incompatible number of variables for merge
   delete default_base_tree;
   default_base_tree= NULL;
   default_base_tree= new tree_type(3, 0, 0, 1);

   delete pda[0];
   pda[0]= NULL;
   delete [] pda;
   pda= new Distribution*[3];
   pda[0]= new Distribution(adistrib);
   pda[1]= pda[0];
   pda[2]= pda[0];

   delete [] itype;
   itype= new int[3];
   itype[0]= INT_VALUE;
   itype[1]= INT_VALUE;
   itype[2]= INT_VALUE;
   v.reset(3, 0);
   pota= new int_trees*[nb_var_trees];
   otrees= new tree_type*[nb_trees];
   for(o= 0; o < nb_var_trees; o++)
   {
      for(t= 0; t < nb_trees; t++)
      {
         tmp_utree->simulation(sdistrib, max_size, max_depth);
         otrees[t]= new tree_type(*default_base_tree);
         otrees[t]->set_structure(*tmp_utree, v);
         otrees[t]->iid_simulation(pda);
         cout << "tree number " << t << ":"
              << "(" << otrees[t]->get_nb_int() << " integral variables,"
              << otrees[t]->get_nb_float() << " floating variables)"<< endl;
         cout << *otrees[t] << endl;
      }
      pota[o]= new Int_trees(nb_trees, itype, otrees);
   }
   cout << "Testing for merge with incompatible number of variables." << endl;
   merged= b->merge(error, nb_var_trees, pota);
   if (error.get_nb_error() > 0)
      cout << error << endl;
   else
   {
      cout << "Failed to detect error(s)." << endl;
      delete merged;
      merged= NULL;
   }

   for(t= 0; t < nb_trees; t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }
   delete [] otrees;
   otrees= NULL;

   delete default_base_tree;
   default_base_tree= NULL;
   default_base_tree= new tree_type(2, 0, 0, 1);
   delete pda[0];
   pda[0]= NULL;
   delete [] pda;
   pda= new Distribution*[2];
   pda[0]= new Distribution(adistrib);
   pda[1]= pda[0];

   delete [] itype;
   itype= new int[2];
   itype[0]= INT_VALUE;
   itype[1]= STATE;
   v.reset(2, 0);
   pota= new int_trees*[nb_var_trees];
   otrees= new tree_type*[nb_trees];
   for(o= 0; o < nb_var_trees; o++)
   {
      for(t= 0; t < nb_trees; t++)
      {
         tmp_utree->simulation(sdistrib, max_size, max_depth);
         otrees[t]= new tree_type(*default_base_tree);
         otrees[t]->set_structure(*tmp_utree, v);
         otrees[t]->iid_simulation(pda);
      }
      pota[o]= new Int_trees(nb_trees, itype, otrees);
   }
   for(t= 0; t < nb_trees; t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }
   delete [] otrees;
   otrees= NULL;

   cout << "Testing for merge with incompatible variable types." << endl;
   merged= b->merge(error, nb_var_trees, pota);
   if (error.get_nb_error() > 0)
      cout << error << endl;
   else
   {
      cout << "Failed to detect error(s)." << endl;
      delete merged;
      merged= NULL;
   }

   // build tree with float variable
   delete default_base_tree;
   default_base_tree= NULL;
   default_base_tree= new tree_type(1, 1, 0, 1);
   rv.reset(1,1);
   // v.reset(1,1);
   ctrees= new rtree_type*[cmerged->get_nb_trees()];
   otrees= new tree_type*[cmerged->get_nb_trees()];
   for(t= 0; t < cmerged->get_nb_trees(); t++)
   {
      otrees[t]= cmerged->get_tree(t);
      tmp_utree= otrees[t]->get_structure();
      ctrees[t]= new rtree_type(*default_base_tree);
      ctrees[t]->set_structure(*tmp_utree, rv);
      // here is the problem
      tie(it, end)= ctrees[t]->vertices();
      while (it < end)
      {
         v= otrees[t]->get(*it);
         rv.Int(0)= v.Int(0);
         rv.Double(0)= v.Int(1)+(double)random() / (double)0x7fffffff;
         ctrees[t]->put(*it++, rv);
      }
   }

   for(t= 0; t < cmerged->get_nb_trees(); t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }
   delete [] otrees;
   otrees= NULL;

   c= new Trees(cmerged->get_nb_trees(), citype, ctrees);
   potc[0]= c;
   d= c->merge_variable(error, 1, potc);
   delete c;
   c= new Trees(*d);
   delete d;
   d= c->shift(error, 1, -1);
   delete c;
   c= d->shift(error, 4, 1.5);
   delete d;

   delete [] potc;
   potc= NULL;

   otrees= new tree_type*[cmerged->get_nb_trees()];
   cout << "Trees with real attributes: " << endl;
   cout << *c;
   for(t= 0; t < c->get_nb_trees(); t++)
   {
      otrees[t]= c->get_tree(t);
      otrees[t]->display(cout, otrees[t]->root());
      cout << endl;
   }

   for(t= 0; t < cmerged->get_nb_trees(); t++)
   {
      delete otrees[t];
      otrees[t]= NULL;
   }
   delete [] otrees;
   otrees= NULL;

   // selection of variables
   delete [] ivariables;
   ivariables= new int[1];
   ivariables[0]= 1;
   cout << "Select variable..." << endl;
   res= c->select_variable(error, 1, ivariables, true);
   cout << error;

   if (res != NULL)
   {
      otrees= new tree_type*[res->get_nb_trees()];
      cout << endl << "Selection of the 1st variable" << endl;
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete [] otrees;
      otrees= NULL;
   }

   delete res;
   res= NULL;

   delete [] ivariables;
   ivariables= new int[2];
   ivariables[0]= 2;
   ivariables[1]= 4;
   res= c->select_variable(error, 2, ivariables, true);
   cout << error;

   if (res != NULL)
   {
      otrees= new tree_type*[res->get_nb_trees()];
      cout << endl << "Selection of the 2nd and 4th variables" << endl;
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete [] otrees;
      otrees= NULL;
      delete res;
      res= NULL;
   }

   delete [] ivariables;
   ivariables= new int[1];
   ivariables[0]= 3;
   res= c->select_variable(error, 1, ivariables, false);
   cout << error;

   if (res != NULL)
   {
      otrees= new tree_type*[res->get_nb_trees()];
      cout << endl << "Discard 3rd variable" << endl;
      res->ascii_write(cout, false);
      for(t= 0; t < res->get_nb_trees(); t++)
      {
         otrees[t]= res->get_tree(t);
         otrees[t]->display(cout, otrees[t]->root());
         cout << endl;
         delete otrees[t];
         otrees[t]= NULL;
      }
      delete [] otrees;
      otrees= NULL;
      delete res;
      res= NULL;
   }

   // test of the StatError messages
   // incompatible number of variables for select_variables

   delete [] ivariables;
   ivariables= new int[4];
   ivariables[0]= 1;
   ivariables[1]= 2;
   ivariables[2]= 3;
   ivariables[3]= 4;
   res= c->select_variable(error, 4, ivariables, false);
   if (error.get_nb_error() > 0)
      cout << error << endl;
   else
   {
      cout << "Failed to detect error(s)." << endl;
      delete res;
      res= NULL;
   }

   delete c;
   c = NULL;

   // building and display of a tree with integral labels
   for(i= 0; i < 10; i++ )
   {
     v.Int(0) = i+1;
     v.Int(1) = 0;
     if ((i == 5) || (i == 7))
        v.Int(1) = 1;
     va[i]= ti.add_vertex(v);
   }

   ti.add_edge(va[0], va[1]);
   ti.add_edge(va[1], va[2]);
   ti.add_edge(va[2], va[3]);
   ti.add_edge(va[2], va[4]);
   ti.add_edge(va[4], va[6]);
   ti.add_edge(va[4], va[5], 1);
   ti.add_edge(va[6], va[7], 1);
   ti.add_edge(va[5], va[8]);
   ti.add_edge(va[5], va[9]);

   assert(ti.is_root(va[0]));

   cout << "A predefined (homemade) tree with edge types... " << endl;
   ti.display(cout, ti.root());

   citrees = new tree_type*[1];
   citrees[0] = &ti;
   citype[1] = INT_VALUE;

   c = new Trees(1, citype, citrees);

   delete [] citrees;
   citrees = NULL;

   // build sequences
   cout << "Build sequences (all possible paths)..." << endl;
   res_seq= c->build_sequences(error, true);
   cout << error;

   if (res_seq != NULL)
   {
      res_seq->ascii_data_write(cout, 'c', false);
      delete res_seq;
      res_seq= NULL;
   }

   cout << "Build sequences (all non redundant paths)..." << endl;
   res_seq= c->build_sequences(error, false);
   cout << error;

   if (res_seq != NULL)
   {
      res_seq->ascii_data_write(cout, 'c', false);
      delete res_seq;
      res_seq= NULL;
   }


   for(t= 0; t < c->get_nb_trees(); t++)
   {
      delete ctrees[t];
      ctrees[t]= NULL;
   }
   delete [] ctrees;
   ctrees= NULL;

   delete b;
   b= NULL;

   delete merged;
   merged= NULL;

   delete c;
   c= NULL;


   for(o= 0; o < nb_var_trees; o++)
   {
      delete pota[o];
      pota[o]= NULL;
   }
   delete [] pota;
   pota= NULL;

   delete pda[0];
   delete [] pda;

   delete [] pota;
   pota= NULL;

   delete tmp_utree;
   delete [] itype;
   delete [] citype;
   delete [] iidentifier;
   delete [] limit;
   delete [] ivariables;
   delete [] symbol_array;

   return 0;

}
