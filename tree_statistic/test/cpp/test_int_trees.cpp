/****************************************************************
 *
 *  Test of the data structure for observed trees as defined
 *  in observed_trees.h
 */

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

using namespace Stat_trees;

int main(void)
{

   typedef Int_trees::tree_type tree_type;
   typedef tree_type::value value;
   const int inf_bound= 1, sup_bound= 4; //n= 5,
   const double probability= 0.6;
   const int ident= UNIFORM;
   const double parameter= D_DEFAULT;
   const int nb_trees= 5;
   const DiscreteParametric distrib(ident, inf_bound, sup_bound, parameter, probability);
   const int max_depth= 3;
   const int max_size= 20;
   int t, *itype;
   value v;
   Unlabelled_typed_edge_tree *tmp_utree;
   tree_type *default_base_tree, tmp_base_tree;
   tree_type  **observed_trees;
   Int_trees *b;
   TreeCharacteristics **tc;
   // Typed_edge_one_int_tree *otrees1;
   DiscreteDistributionData *mhisto, *ihsize, *ihnb_children;
   StatError error;


   v.reset(1,0);
   v.Int(0)= 1;

   cout << "Tree attributes : " << v << endl;
   // We check that the basic structures work

   itype= new int[1];
   observed_trees= new tree_type*[nb_trees];
 //    observed_trees= new Int_fl_trees[nb_trees];
 //    not equivalent at all !

   cout << "Default tree : " << endl;
   default_base_tree= new tree_type(1, 0, 0, 1);
   // (nb_integral, nb_float, root, nb_vertices)
   default_base_tree->add_vertex(v);
   // necessary : default_base_tree was empty at previous stage
   default_base_tree->display(cout, default_base_tree->root());

   tmp_utree= new Unlabelled_typed_edge_tree;
   cout << "Default structure : " << endl;
   tmp_utree->display(cout, tmp_utree->root());

   itype[0]= INT_VALUE;
   cout << "itype[0] : " << itype[0] << endl;
   observed_trees[0]= new tree_type(*default_base_tree);
   cout << "*(observed_trees[0]) : " << endl;
   observed_trees[0]->display(cout, observed_trees[0]->root());

   delete observed_trees[0];
   // Random assignement of the structures

   //cout << "Simulated random structures : " << endl;
   for(t= 0; t < nb_trees; t++)
   {
      tmp_utree->simulation(distrib, max_size, max_depth);
      //tmp_utree->display(cout, tmp_utree->root());
      observed_trees[t]= new tree_type(*default_base_tree);
      observed_trees[t]->set_structure(*tmp_utree,v);
      //cout << "(size " << observed_trees[t]->get_size() << ", depth "
      //     << observed_trees[t]->get_depth() << ")" << endl;
      //cout << endl;
   }

   b= new Int_trees(nb_trees, itype, observed_trees);

   b->iid_simulation(distrib);

   // Random assignement of the labels

   cout << "Observed_trees has " << b->get_nb_trees() << " trees..." << endl
        << "... of maximal size " << b->get_max_size() << " and maximal depth "
        << b->get_max_depth() << endl;

   // number of integral and float variables

   cout << "Number of integral variables == " << b->get_nb_int() << endl;
   cout << "Number of float variables == " << b->get_nb_float() << endl;

   // min and max values
   cout << "Min value == " << b->get_min_int_value(0) << endl;
   cout << "Max value == " << b->get_max_int_value(0) << endl;

   for(t= 0; t < nb_trees; t++)
       (*observed_trees[t])= *(b->get_tree(t));

   cout << "Simulated random labels : " << endl;
   for(t= 0; t < nb_trees; t++)
   {
       observed_trees[t]->display(cout, observed_trees[t]->root());
       cout << endl;
   }

   // test of the ascii_write procedure
   cout << "Ascii print of Observed trees (summary): " << endl;
   cout << *b; // b->ascii_write(cout, true);
   cout << endl;

   cout << "Line write of Observed trees : " << endl;
   b->line_write(cout);
   cout << endl << endl;

   cout << "Ascii print of Observed trees (detailed): " << endl;
   b->ascii_write(cout, true);
   cout << endl;


   cout << "Histograms for the tree size : " << endl;
   b->ascii_write_size_frequency_distribution(cout, 1, 0);
   cout << endl;
   cout << "Histograms for the number of children  : " << endl;
   b->ascii_write_nb_children_frequency_distribution(cout, 1, 0);
   cout << endl;

   ihsize= b->extract_size();
   ihnb_children= b->extract_nb_children();

   tc= new TreeCharacteristics*[b->get_nb_int()];
   if (b->get_characteristics(0) != NULL)
   {
      tc[0]= b->get_characteristics(0);

      // delete b;

      cout << "Histogram for the first occurrence of each value : " << endl;
      tc[0]->ascii_write_first_occurrence_root(cout, 1, 0);

      delete tc[0];
   }
   else
   {
      tc[0]= NULL;
      cout << "Get_characteristics() returned a NULL pointer." << endl;
   }

   // Constructor of Int_trees using histograms

   delete [] tc;
   delete b;

   b= new Int_trees(1, *ihsize, *ihnb_children);
   cout << "Constructor of Int_trees using histograms." << endl;
   cout << "Simulated structures : " << endl;
   for(t= 0; t < nb_trees; t++)
       (*observed_trees[t])= *(b->get_tree(t));
   for(t= 0; t < nb_trees; t++)
   {
       observed_trees[t]->display(cout, observed_trees[t]->root());
       cout << endl;
   }

   cout << "Histograms for the tree size : " << endl;
   b->ascii_write_size_frequency_distribution(cout, 1, 0);
   cout << endl;
   cout << "Histograms for the number of children  : " << endl;
   b->ascii_write_nb_children_frequency_distribution(cout, 1, 0);
   cout << endl;

   delete b;

   // Same principe for 2-dimensionnal trees

   // We check that the basic structure work

   itype= new int[2];
   itype[0]= INT_VALUE;
   itype[1]= INT_VALUE;
   // Random assignement of the structures

   v.reset(2,0);
   cout << endl << "Tree attributes : " << v << endl;
   // We check that the basic structure work

   cout << "Default tree : " << endl;
   default_base_tree= new tree_type(2, 0, 0, 1);
   default_base_tree->add_vertex(v);
   // necessary : default_base_tree was empty at previous stage
   default_base_tree->display(cout,default_base_tree->root());

   for(t= 0; t < nb_trees; t++)
   {
      tmp_utree->simulation(distrib, max_size, max_depth);
      //tmp_utree->display(cout, tmp_utree->root());
      *(observed_trees[t])= *default_base_tree;
      observed_trees[t]->set_structure(*tmp_utree,v);
      //cout << endl;
   }

   b= new Int_trees(nb_trees, itype, observed_trees);

   cout << "Observed_trees has " << b->get_nb_trees() << " trees." << endl;

   b->iid_simulation(distrib);

   // Random assignement of the labels
   // number of integral and float variables

   cout << "Number of integral variables == " << b->get_nb_int() << endl;
   cout << "Number of float variables == " << b->get_nb_float() << endl;

   // min and max values
   cout << "Min value == " << b->get_min_int_value(0) << " (1st variable)" << endl;
   cout << "Max value == " << b->get_max_int_value(0) << " (1st variable)" << endl;

   cout << "Min value == " << b->get_min_int_value(1) << " (2nd variable)" << endl;
   cout << "Max value == " << b->get_max_int_value(1) << " (2nd variable)" << endl;

   cout << "Simulated random labels : " << endl;

   for(t= 0; t < nb_trees; t++)
   {
       *(observed_trees[t])= *(b->get_tree(t));
       observed_trees[t]->display(cout, observed_trees[t]->root());
       cout << "(size " << observed_trees[t]->get_size() << ", depth "
            << observed_trees[t]->get_depth() << ")" << endl;
       cout << endl;
   }

   tc= new TreeCharacteristics*[b->get_nb_int()];
   if (b->get_characteristics(0) != NULL)
   {
      tc[0]= b->get_characteristics(1);

      cout << "Histogram for the marginal distribution (2nd variable) : " << endl;
      tc[0]->ascii_write_marginal_distribution(cout, 1, 0);

      delete tc[0];
   }
   else
   {
      tc[0]= NULL;
      cout << "Get_characteristics() returned a NULL pointer." << endl;
   }



   // test of the copy constructor and assignement operator
   // of Int_trees

   Int_trees *pt_b= new Int_trees(*b), cp_b= *pt_b;
   delete b;

   cout << endl << "Test of the copy constructor of 'Int_trees' " << endl;
   cout << "Number of integral variables == " << pt_b->get_nb_int() << endl;
   cout << "Number of float variables == " << pt_b->get_nb_float() << endl;

   // min and max values
   cout << "Min value == " << pt_b->get_min_int_value(0) << endl;
   cout << "Max value == " << pt_b->get_max_int_value(0) << endl;

   cout << "Number of observed trees : " << pt_b->get_nb_trees() << endl;
   cout << "Last observed tree : " << endl;
   (pt_b->get_tree(pt_b->get_nb_trees()-1))->display(cout, 0);

   cout << "Histogram for marginal distribution " << endl;
   (pt_b->get_characteristics(0))->ascii_write_marginal_distribution(cout, 1, 0);
   delete pt_b;

   cout << endl << "Test of the assignement operator of 'Int_trees' " << endl;
   cout << "Number of integral variables == " << cp_b.get_nb_int() << endl;
   cout << "Number of float variables == " << cp_b.get_nb_float() << endl;

   // min and max values
   cout << "Min value == " << cp_b.get_min_int_value(0) << endl;
   cout << "Max value == " << cp_b.get_max_int_value(0) << endl;

   cout << "Number of observed trees : " << cp_b.get_nb_trees() << endl;
   cout << "Last observed tree : " << endl;
   (cp_b.get_tree(cp_b.get_nb_trees()-1))->display(cout, 0);

   cout << "Histogram for marginal distribution (1st variable)" << endl;
   (cp_b.get_characteristics(0))->ascii_write_marginal_distribution(cout, 1, 0);

   mhisto= cp_b.extract(error, 1);
   cout << "Distribution_data for marginal distribution (1st variable)" << endl;
   cout << *mhisto;

   for(t= 0; t < nb_trees; t++)
      delete observed_trees[t];

   delete [] tc;
   delete [] itype;

   delete [] observed_trees;
   delete default_base_tree;
   delete tmp_utree;
   delete mhisto;

   return 0;

}
