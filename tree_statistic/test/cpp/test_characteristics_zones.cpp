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

  typedef Int_trees::tree_type  tree_type;
  typedef tree_type::value value;
  typedef tree_type::vertex_array vertex_array;

  const int n= 10, inf_bound= 0, sup_bound= 3;
  const double probability= 0.6;
  const int ident= UNIFORM;
  const double parameter= D_DEFAULT;
  int nb_trees= n;
  const DiscreteParametric cdistrib0(ident,inf_bound,sup_bound,parameter,probability);
  const DiscreteParametric cdistribl(ident,2,5,parameter,probability);
  int cmax_depth= 20;
  int cmax_size= 20;
  int t, *itype;
  int _nb_integral;
  unsigned int _max_value, _min_value;
  unsigned int _max_depth, _max_size, i;
  value v;
  Unlabelled_typed_edge_tree *tmp_utree;
  tree_type *default_base_tree, tmp_base_tree;
  tree_type  **observed_trees;
  Int_trees *b;
  Typed_edge_one_int_tree **otrees1;
  TreeCharacteristics **tc;
  Distribution **distrib;
  FrequencyDistribution h;
  vertex_array va;
  generic_visitor<tree_type> visitor;

  v.reset(2,0);
  v.Int(0)= 1;
  _nb_integral= 1;

  distrib= new Distribution*[2];
  distrib[0]= new Distribution(cdistrib0);
  distrib[1]= new Distribution(cdistribl);
  // distribution used for the label simulation

  itype= new int[2];
  observed_trees= new tree_type*[nb_trees];
  otrees1= new Typed_edge_one_int_tree*[nb_trees];

  default_base_tree= new tree_type(2,0,0,1);
  default_base_tree->add_vertex(v);
  // necessary : default_base_tree is empty at this stage

  tmp_utree= new Unlabelled_typed_edge_tree;

  itype[0]= INT_VALUE;
  itype[1]= INT_VALUE;
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
     va= visitor.get_breadthorder(*(observed_trees[t]));
     for(i= 0; i < observed_trees[t]->get_size(); i++)
     {
        v.Int(0)= (observed_trees[t]->get(va[i])).Int(0);
        v.Int(1)= va[i];
        observed_trees[t]->put(va[i], v);
     }
     i= (observed_trees[t]->get_min_value())->Int(0);
     _min_value= min(_min_value, i);
     i= (observed_trees[t]->get_max_value())->Int(0);
     _max_value= max(_max_value, i);
     _max_depth= max(_max_depth, observed_trees[t]->get_depth());
     _max_size= max(_max_size, observed_trees[t]->get_size());
  }

  // build the histogramms
  cout << "_min_value == " << _min_value << ", _max_value == " << _max_value
       << " _max_depth == " << _max_depth << ", _max_size == " << _max_size << endl;

  tc= new TreeCharacteristics*[_nb_integral];

  for (t= 0; t < nb_trees; t++)
      otrees1[t]= observed_trees[t]->select_int_variable(0);


  tc[0]= new TreeCharacteristics(_min_value,
                                  _max_value,
                                  _max_size,
                                  _max_depth,
                                  nb_trees,
                                  otrees1,
                                  0);

  cout << "Histogram for the size of the homogeneous zones : " << endl;

  tc[0]->ascii_write_sojourn_size(cout, 1, 0);

  cout << "Histogram for the number of homogeneous zones : " << endl;

  tc[0]->ascii_write_nb_zones(cout, 1, 0);

  b= new Int_trees(nb_trees, itype, observed_trees);

  cout << "Observed_trees has " << b->get_nb_trees() << " trees of maximal size "
       << b->get_max_size() << " and maximal depth " << b->get_max_depth() << endl;
  cout << "Simulated random labels : " << endl;

  for(t= 0; t < nb_trees; t++)
  {
      (*observed_trees[t])= *(b->get_tree(t));
      observed_trees[t]->display(cout, observed_trees[t]->root());
      cout << endl;
  }

  // number of integral and float variables

  cout << "Number of integral variables == " << b->get_nb_int() << endl;
  cout << "Number of float variables == " << b->get_nb_float() << endl;

  // min and max values
  cout << "Min value == " << b->get_min_int_value(0) << endl;
  cout << "Max value == " << b->get_max_int_value(0) << endl;

  delete b;

  for(t= 0; t < nb_trees; t++)
     delete observed_trees[t];

  delete [] observed_trees;
  delete [] otrees1;
  delete tc[0];
  delete [] tc;
  delete distrib[0];
  delete [] distrib;
  delete [] itype;
  delete default_base_tree;
  delete tmp_utree;

  return 0;

}
