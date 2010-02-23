/****************************************************************
 *
 *  Test of the int_fl_tree methods as defined in observed_trees.h
 *  and observed_int_trees.h
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

  typedef tree_type::vertex_descriptor vertex;
  typedef tree_type::value value;
  typedef tree_type::key key;
  typedef tree_type::children_iterator children_iterator;
  typedef Unlabelled_typed_edge_tree::vertex_descriptor char_vertex;
  typedef Unlabelled_typed_edge_tree::children_iterator char_children_iterator;

  const int n= 5, inf_bound= 0, sup_bound= 3, _nb_trees= 2;
  const double probability= 0.6;
  const int ident= UNIFORM;
  const double parameter= D_DEFAULT;
  // int max_depth= 3;
  // int nv= 15;

  const DiscreteParametric cdistrib0(ident,inf_bound,sup_bound,parameter,probability);
  const DiscreteParametric cdistrib1(ident,1,5,parameter,probability);
  Distribution **distrib;

  tree_type *pt, t, **trees;
  Typed_edge_one_int_tree o, *otrees1, *tmp_otree1;
  Unlabelled_typed_edge_tree utree, *tmp_utree;
  // default constructor;
  DiscreteParametric var_distrib(BINOMIAL,inf_bound,sup_bound,parameter,probability);
  value v, val;
  vertex *vt= new vertex[n], child;
  int i, var;

  pt= new tree_type();
  distrib= new Distribution*[2];

  distrib[0]= new DiscreteParametric(cdistrib0);
  distrib[1]= new DiscreteParametric(cdistrib1);

  // Default Int_fl_tree

  pt->add_vertex(v);
  cout << "Default Int_fl_tree : " << endl;
  pt->display(cout, pt->root());

  delete pt;

  // One int Int_fl_tree

  pt=new tree_type(1,0,0,5);
  v.reset(1,0);
  pt->add_vertex(v);
  cout << "Int_fl_tree with " << pt->get_nb_int() << " int and " <<
       pt->get_nb_float() << " double : " << endl;
  pt->display(cout, pt->root());
  v.Int(0)= 1;
  pt->put(0,v);
  vt[0]= 0;

  for(i= 1; i < n; i++ )
  {
     v.Int(0)= i+1;
     vt[i]= pt->add_vertex(v);
  }

  // edges
  pt->add_edge(vt[0], vt[1]);
  pt->add_edge(vt[0], vt[2]);
  pt->add_edge(vt[1], vt[3]);
  pt->add_edge(vt[1], vt[4]);

  assert(pt->is_root(vt[0]));

  cout << "Tree set to : " << endl;
  pt->display(cout, pt->root());

  // Copy constructeur
  tree_type cpt(*pt);

  cout << "A copy of this tree (constructor) : " << endl;
  cpt.display(cout,cpt.root());

  t= cpt;

  cout << "A copy of this tree (assignement) : " << endl;
  t.display(cout,t.root());

  delete pt;

  // Two int one float Int_fl_tree

  pt=new tree_type(2,1,0,5);
  v.reset(2,1);
  pt->add_vertex(v);
  cout << endl;
  cout << "Int_fl_tree with " << pt->get_nb_int() << " int and " <<
       pt->get_nb_float() << " double : " << endl;
  pt->display(cout, pt->root());
  v.Int(0)= 1;
  v.Int(1)= 2;
  v.Double(0)= v.Int(0)/2.0;
  pt->put(0,v);
  vt[0]= 0;

  for(i= 1; i < n; i++ )
  {
     v.Int(0)= i+1;
     v.Int(1)= i+2;
     v.Double(0)= v.Int(0)/2.0;
     vt[i]= pt->add_vertex(v);
  }

  // edges
  pt->add_edge(vt[0], vt[1]);
  pt->add_edge(vt[0], vt[2]);
  pt->add_edge(vt[1], vt[3]);
  pt->add_edge(vt[1], vt[4]);

  assert(pt->is_root(vt[0]));

  cout << "Tree set to : " << endl;
  pt->display(cout, pt->root());

  // Copy constructeur
  tree_type c2pt(*pt);

  cout << "A copy of this tree (constructor) : " << endl;
  cout << "with " << c2pt.get_nb_int() << " int and " <<
       c2pt.get_nb_float() << " double : " << endl;

  c2pt.display(cout,c2pt.root());

  t= c2pt;

  cout << "A copy of this tree (assignement) : " << endl;
  cout << "with " << t.get_nb_int() << " int and " <<
       t.get_nb_float() << " double : " << endl;
  t.display(cout,t.root());

  // variable selection

  o= *(t.select_int_variable(1));

  cout << "Selection of the 2nd variable of this tree : " << endl;
  o.display(cout,o.root());

  trees= new tree_type*;
  *trees= new tree_type[_nb_trees];
  otrees1= new Typed_edge_one_int_tree[_nb_trees];

  (*trees)[0]= t;
  (*trees)[1]= c2pt;

  cout << _nb_trees << " tree(s) with " << t.get_nb_int() << " int and "
       << t.get_nb_float() << " float " << endl;
  for (var= 0; var < t.get_nb_int(); var++)
  {
      for (i= 0; i < _nb_trees; i++)
      {
          tmp_utree= (*trees)[i].get_structure();
          tmp_otree1= (*trees)[i].select_int_variable(var);
          otrees1[i]= *tmp_otree1;
          cout << "Selection of the " << var+1 << "th variable of the tree : " << endl;
          otrees1[i].display(cout, otrees1[i].root());
      }
  }

  utree.simulation(*distrib[0],20,3);
# ifdef DEBUG
  cout << "Tree used for simulation" << endl;
  utree.display(cout, utree.root());
# endif

  // random simulation using the unlabelled-tree constructors
  pt= new tree_type(utree, v);
  cout << endl;
  cout << "Random-structure tree with " << pt->get_nb_int() << " int and " <<
       pt->get_nb_float() << " double : " << endl;
  pt->display(cout, pt->root());

  delete pt;

  pt= new tree_type(2, 0, utree);

  cout << endl;
  cout << "Random-structure tree with " << pt->get_nb_int() << " int and " <<
       pt->get_nb_float() << " double : " << endl;
  pt->display(cout, pt->root());

  cout << "The above tree has " << pt->get_size() << " vertices and depth "
       << pt->get_depth() << endl;

  cout << "Random affectation of the labels " << endl;
  pt->iid_simulation(distrib);
  pt->display(cout, pt->root());

  cout << "Its maximal values are : " << *(pt->get_max_value())
       << " and its minimal values are : " << *(pt->get_min_value()) << endl;

  t= *pt;
  child= 2;
  val= t.get(child);
  delete pt;
  pt= NULL;
  pt= select_subtree(t, child);

  if (pt != NULL)
  {
     cout << "Extracting subtree rooted at node " << val
          << "from this tree " << endl;
     pt->display(cout, pt->root());
     delete pt;
     pt= NULL;
  }

  pt= select_subtree(t, child, false);

  if (pt != NULL)
  {
     cout << "Pruning this tree at node " << val << endl;
     pt->display(cout, pt->root());
     delete pt;
     pt= NULL;
  }

  delete [] vt;
  delete distrib[0];
  delete distrib[1];
  delete [] distrib;
  delete tmp_utree;
  delete tmp_otree1;
  delete [] otrees1;

  return 0;

}
