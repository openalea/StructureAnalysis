/****************************************************************
 *
 *  Test of the data structure for observed trees as defined
 *  in observed_trees.h
 */

#include "TREE/tree_simple.h"
#include "TREE/tree_traits.h"
#include "TREE/basic_visitors.h"
#include "stat_tools.h"
#include "generic_typed_edge_tree.h"
#include "int_fl_containers.h"
#include "curves.h"
#include "markovian.h"
#include "sequences.h"
#include "typed_edge_trees.h"
#include "int_trees.h"

using namespace Stat_trees;

int main(void)
{

  typedef Int_trees::tree_type  tree_type;
  typedef tree_type::value value;
  const int n= 10, inf_bound= 0, sup_bound= 3;
  const double probability= 0.6;
  const int ident= UNIFORM;
  const double parameter= D_DEFAULT;
  int nb_trees= n;
  const Parametric cdistrib0(ident,inf_bound,sup_bound,parameter,probability);
  const Parametric cdistribl(ident,2,5,parameter,probability);
  int cmax_depth= 20;
  int cmax_size= 20;
  int t, *itype, i;
  int _nb_integral, _max_value, _min_value;
  unsigned int _max_depth, _max_size;
  value v;
  Unlabelled_typed_edge_tree *tmp_utree;
  tree_type *default_base_tree, tmp_base_tree;
  tree_type  **observed_trees;
  // Int_trees *b;
  Typed_edge_one_int_tree *otrees1;
  Tree_characteristics **tc;
  Distribution **distrib;
  Histogram h;
  v.reset(1,0);
  v.Int(0)= 1;
  _nb_integral= 1;

  distrib= new Distribution*[0];
  distrib[0]= new Distribution(cdistribl);
  // distribution for the label simulation

  itype= new int[0];
  observed_trees= new tree_type*[nb_trees];
  otrees1= new Typed_edge_one_int_tree[nb_trees];

  default_base_tree= new tree_type(1,0,0,1);
  default_base_tree->add_vertex(v);
  // necessary : default_base_tree is empty at this stage

  tmp_utree= new Unlabelled_typed_edge_tree;

  itype[0]= INT_VALUE;
  // Random assignement of the structures

  _min_value= INT_MAX;
  _max_value= 0;
  _max_depth= 0;
  _max_size= 1;

  for(t= 0; t < nb_trees; t++)
  {
     tmp_utree->simulation(cdistrib0, cmax_size, cmax_depth);
     observed_trees[t]= new tree_type(*default_base_tree);
     observed_trees[t]->set_structure(*tmp_utree,v);
     observed_trees[t]->iid_simulation(distrib);
     i= (observed_trees[t]->get_min_value())->Int(0);
     _min_value= min(_min_value, i);
     i= (observed_trees[t]->get_max_value())->Int(0);
     _max_value= max(_max_value, i);
     _max_depth= max(_max_depth, observed_trees[t]->get_depth());
     _max_size= max(_max_size, observed_trees[t]->get_size());
  }

  // building of the histogramms
  cout << "_min_value == " << _min_value << ", _max_value == " << _max_value
       << " _max_depth == " << _max_depth << ", _max_size == " << _max_size << endl;

  tc= new Tree_characteristics*[_nb_integral];

  for (t= 0; t < nb_trees; t++)
      otrees1[t]= *(observed_trees[t]->select_int_variable(0));


  tc[0]= new Tree_characteristics(_min_value,
                                  _max_value,
                                  _max_size,
                                  _max_depth,
                                  nb_trees,
                                  otrees1,
                                  0);


  // copy of the Tree_characteristics objects;

  cout << "Checking the assignement operator of 'Tree_characteristics' : " << endl;
  Tree_characteristics *pt_copy_tc= new Tree_characteristics(*(tc[0])), copy_tc= *(tc[0]);

  cout << "Histogram for the marginal distribution is : " << endl;
  tc[0]->ascii_write_marginal(cout, 1, 0);

  //delete tc[0];
  //delete [] tc;

  //cout << "Histogram of the copy (marginal) : " << endl;
  //copy_tc.ascii_write_marginal(cout, 1, 0);

  cout << "Histogram for the number of occurrences for each value : " << endl;
  tc[0]->ascii_write_nb_occurrences(cout, 1, 0);

  //cout << "Histogram of the copy (occurrences) : " << endl;
  //copy_tc.ascii_write_nb_occurrences(cout, 1, 0);

  cout << "Histogram for the number of zones for each value : " << endl;
  tc[0]->ascii_write_nb_zones(cout, 1, 0);

  //cout << "Histogram of the copy (zones) : " << endl;
  //copy_tc.ascii_write_nb_zones(cout, 1, 0);

  //cout << "Checking the copy constructor of 'Tree_characteristics' : " << endl;

  //cout << "Histogram of the copy (marginal) : " << endl;
  //pt_copy_tc->ascii_write_marginal(cout, 1, 0);

  //cout << "Histogram of the copy (occurrences) : " << endl;
  //pt_copy_tc->ascii_write_nb_occurrences(cout, 1, 0);

  //cout << "Histogram of the copy (zones) : " << endl;
  //pt_copy_tc->ascii_write_nb_zones(cout, 1, 0);

  delete pt_copy_tc;

  for(t= 0; t < nb_trees; t++)
     delete observed_trees[t];

  delete [] observed_trees;
  delete [] otrees1;
  delete distrib[0];
  delete [] distrib;
  delete default_base_tree;
  delete tmp_utree;

  return 0;

}
