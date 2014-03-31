/****************************************************************
 *
 *  Test of the generic_tree methods as defined in generic_tree.h
 */

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"
#include "stat_tool/stat_tools.h"
#include "tree_statistic/generic_typed_edge_tree.h"
#include <iostream>
#include <string>
#include <deque>
#include <cstdlib>

using namespace Stat_trees; 
using namespace Tree_tie;

typedef Generic_tree<int> tree_int;

int main(void)
{
   const int inf_bound= 0, sup_bound= 3;
   const double probability= 0.6;
   const int ident= UNIFORM;
   const double parameter= D_DEFAULT;
   int max_depth= 3;
   int nv= 15; //n= 5,
   std::vector<int> depth;
   unsigned int i, n= 5;

   // const DiscreteParametric distrib(nb_value+3,ident);
   const DiscreteParametric distrib(ident,inf_bound,sup_bound,parameter,probability);


   tree_int t;
   // default constructor;
   typedef tree_int::vertex_descriptor vertex;
   typedef tree_int::key key;
   typedef tree_int::children_iterator children_iterator;
   typedef tree_int::value value;
   typedef generic_visitor<tree_int>::vertex_array vertex_array;
   typedef generic_visitor<tree_int>::vertex_deque vertex_deque;
   typedef Unlabelled_typed_edge_tree::vertex_descriptor char_vertex;
   typedef Unlabelled_typed_edge_tree::children_iterator char_children_iterator;
   DiscreteParametric var_distrib(BINOMIAL,inf_bound,sup_bound,parameter,probability);
   vertex *v= new vertex[n];


   // building and display of a tree with integral labels
   for(i= 0; i < n; i++ )
     v[i]= t.add_vertex(i+1);

   // edges
   t.add_edge(v[0], v[1]);
   t.add_edge(v[0], v[2]);
   t.add_edge(v[1], v[3]);
   t.add_edge(v[1], v[4]);

   assert(t.is_root(v[0]));

   tree_int::vertex_iterator it, end;
   // key root= t.root();
   children_iterator ci, cend;
   tie(ci, cend) = t.children(t.root());
   vertex child= *ci++;

   tie(it, end)= t.vertices();
   while (it < end)
      t.put(*it++, *it+1);

   cout << "A predefined (homemade) tree ... " << endl;
   t.display(cout, t.root());

   cout << "... of depth " << t.get_depth() << " and size " <<
        t.get_size() << "." << endl;

   for(i= 0; i < n; i++)
      cout << "Number of children of node " << i+1 << ": " << t.get_nb_children(i)
           << "; depth: " << t.get_depth(i) << endl;
   cout << endl;

   cout << "Subtree from this tree rooted at node " << child+1 << endl;

   t.display(cout, child);

   generic_visitor<tree_int> *visitor= new generic_visitor<tree_int>;
   traverse_tree(t.root(),t,*visitor);

   visitor->print_preorder(t);
   visitor->print_inorder(t);
   visitor->print_postorder(t);

   cout << "Breath first tree traversal: " << endl;
   vertex_array va= visitor->get_breadthorder(t);
   for(i= 0; i < va.size(); i++)
      cout << va[i]+1 << endl;

   cout << endl << "Leaves last tree traversal : " << endl;
   va= visitor->get_leavesfirstorder(t, depth);
   for(i= 0; i < va.size(); i++)
      cout << va[i]+1 << endl;
   cout << endl;

   cout << "Distance of each node to nearest leaf : " << endl;
   for(i= 0; i < va.size(); i++)
      cout << depth[i] << endl;
   cout << endl;

   cout << "Ancestors of vertex 5: " << endl;
   vertex_deque vd= visitor->get_vertex_ancestors(t, 4);
   for(i= 0; i < vd.size(); i++)
      cout << vd[i]+1 << endl;
   cout << endl;

   // extracting a subtree
   tree_int *ptcopy_t= NULL;
   ptcopy_t= select_subtree(t, child);

   if (ptcopy_t != NULL)
   {
      cout << "Extracting subtree from this tree rooted at node " << child+1 << endl;
      ptcopy_t->display(cout, ptcopy_t->root());
      delete ptcopy_t;
   }
   ptcopy_t= NULL;

   // extracting a trivial subtree
   // vertices are renumbered from 0
   // but labels are copied
   ptcopy_t= select_subtree(t, 3);

   if (ptcopy_t != NULL)
   {
      cout << "Extracting subtree from this tree rooted at node " << 4 << endl;
      ptcopy_t->display(cout, ptcopy_t->root());
      cout << "New id of root node: " << ptcopy_t->root() << endl;
      visitor= new generic_visitor<tree_int>;
      traverse_tree(ptcopy_t->root(), *ptcopy_t, *visitor);

      cout << "Breath first tree traversal starting at root vertex: " << endl;
      va= visitor->get_breadthorder(*ptcopy_t, ptcopy_t->root());
      for(i= 0; i < va.size(); i++)
         cout << va[i]+1 << endl;
      delete ptcopy_t;
      delete visitor;
      visitor= NULL;
   }
   ptcopy_t= NULL;

   ptcopy_t= select_subtree(t, child, false);

   if (ptcopy_t != NULL)
   {
      cout << "Pruning this tree at node " << child+1 << endl;
      ptcopy_t->display(cout, ptcopy_t->root());
      delete ptcopy_t;
   }
   ptcopy_t= NULL;

   delete visitor;
   visitor= new generic_visitor<tree_int>;
   traverse_tree(child, t, *visitor);

   cout << "Preorder tree traversal starting at vertex "
        << child+1 << ": " << endl;
   visitor->print_preorder(t);

   cout << "Breath first tree traversal starting at vertex "
        << child+1 << ": " << endl;
   va= visitor->get_breadthorder(t, child);
   for(i= 0; i < va.size(); i++)
      cout << va[i]+1 << endl;

   // construction of a tree by copy
   ptcopy_t= new tree_int(t);
   tie(it,end)= ptcopy_t->vertices();
   while( it != end )
            ptcopy_t->put(*it++,*it+1);

   cout << "A copy of the whole tree is (constr) " << endl;
   ptcopy_t->display(cout, ptcopy_t->root());
   cout << "It has depth " << ptcopy_t->get_depth() << " and size " <<
        ptcopy_t->get_size() << "." << endl;


   tree_int copy_t;
   copy_t= *ptcopy_t;
   delete ptcopy_t;
   tie(it,end)= copy_t.vertices();
   while( it != end )
            copy_t.put(*it++,*it+1);

   cout << "Another copy : " << endl;
   copy_t.display(cout, copy_t.root());
   cout << "... which has depth " << copy_t.get_depth() << " and size " <<
        copy_t.get_size() << "." << endl;


   ptcopy_t= new tree_int(6, 10);

   delete [] v;
   delete visitor;
   n= 10;
   v= new vertex[n];
   visitor= new generic_visitor<tree_int>;

   // building and display of a tree with integral labels
   for(i= 0; i < n; i++ )
     v[i]= ptcopy_t->add_vertex(i+1);

   // edges
   ptcopy_t->add_edge(v[6], v[5]);
   ptcopy_t->add_edge(v[6], v[2]);
   ptcopy_t->add_edge(v[6], v[4]);
   ptcopy_t->add_edge(v[5], v[8]);
   ptcopy_t->add_edge(v[5], v[1]);
   ptcopy_t->add_edge(v[5], v[0]);
   ptcopy_t->add_edge(v[4], v[3]);
   ptcopy_t->add_edge(v[4], v[7]);
   ptcopy_t->add_edge(v[0], v[9]);

   tie(it, end)= ptcopy_t->vertices();
   while (it < end)
      ptcopy_t->put(*it++, *it);

   assert(ptcopy_t->is_root(v[6]));

   cout << "A second predefined tree ... " << endl;
   ptcopy_t->display(cout, ptcopy_t->root());

   cout << "... of depth " << ptcopy_t->get_depth() << " and size " <<
        ptcopy_t->get_size() << "." << endl;

   cout << "Subtree from this tree rooted at node " << 4 << endl;

   ptcopy_t->display(cout, 4);
   t= *ptcopy_t;
   delete ptcopy_t;
   value val;

   child= 5;
   val= t.get(child);
   ptcopy_t= select_subtree(t, child);

   if (ptcopy_t != NULL)
   {
      cout << "Extracting subtree from this tree rooted at node " << val << endl;
      ptcopy_t->display(cout, ptcopy_t->root());
      delete ptcopy_t;
   }
   ptcopy_t= NULL;
   ptcopy_t= select_subtree(t, child, false);

   if (ptcopy_t != NULL)
   {
      cout << "Pruning this tree at node " << val << endl;
      ptcopy_t->display(cout, ptcopy_t->root());
      delete ptcopy_t;
   }

   ptcopy_t= new tree_int(t);
   traverse_tree(ptcopy_t->root(),*ptcopy_t,*visitor);

   visitor->print_preorder(*ptcopy_t);
   visitor->print_inorder(*ptcopy_t);
   visitor->print_postorder(*ptcopy_t);

   cout << "Breath first tree traversal: " << endl;
   va= visitor->get_breadthorder(*ptcopy_t);
   for(i= 0; i < va.size(); i++)
      cout << va[i] << endl;

   // for(i= 0; i < n; i++)
   //    cout << "Number of children of node " << va[i] << " : "
   //         << ptcopy_t->get_nb_children(va[i]) << endl;
   // cout << endl;

   cout << endl << "Leaves-last tree traversal: " << endl;
   va= visitor->get_leavesfirstorder(*ptcopy_t, depth);
   for(i= 0; i < va.size(); i++)
      cout << va[i] << endl;
   cout << endl;

   cout << "Distance of each node to nearest leaf : " << endl;
   for(i= 0; i < va.size(); i++)
      cout << depth[i] << endl;
   cout << endl;

   delete [] v;
   delete ptcopy_t;
   n= 5;
   // building and display of an unlabelled tree
   // i.e. a tree structure

   cout << "Tree structures have label : " << C_DEFAULT_CHAR << endl << endl;

   // direct building of a unlabelled tree
   Unlabelled_typed_edge_tree utree; // (0, n);
   char_vertex *char_v= new char_vertex[n];
   char_children_iterator char_ci, char_cend;

   // root= utree.add_vertex(C_DEFAULT_CHAR);
   // utree.root()= root;
   char_v[0]= utree.root();

   for(i= 1; i < n; i++ )
        char_v[i]= utree.add_vertex(C_DEFAULT_CHAR);

   assert(utree.is_root(char_v[0]));

   utree.add_edge(char_v[0], char_v[1]);
   utree.add_edge(char_v[0], char_v[2]);
   utree.add_edge(char_v[1], char_v[3]);
   utree.add_edge(char_v[1], char_v[4]);

   tie(it,end)= utree.vertices();
   cout << "Vertices are : " << *it+1 << " ... " << *(end-1)+1 << endl;
   while( it != end )
            utree.put(*it++,C_DEFAULT_CHAR);
   cout << "Root node is : " << utree.root()+1 << endl;

   tie(char_ci, char_cend)= utree.children(utree.root());
 //  assert(utree.is_leaf(*char_ci+1));

 //  cout << "Root node has children : " << *char_ci << endl;

                     // " << utree.root() <<
   cout << "... Above tree has following structure (built manually) " << endl;
   utree.display(cout, utree.root());
   cout << "It has depth " << utree.get_depth() << " and size " <<
        utree.get_size() << "." << endl;

   // creation of a tree structure from a tree structure
   Unlabelled_typed_edge_tree copy_utree(utree);

   cout << "... whose copy is : " << endl;
   copy_utree.display(cout, copy_utree.root());
   cout << "(" << copy_utree.get_size() << " vertices and depth "
        << copy_utree.get_depth() << ")" << endl;

   // creation of a tree structure from a generic tree
   Unlabelled_typed_edge_tree *pt_utree= new Unlabelled_typed_edge_tree(copy_t);

   cout << "This can also be built from the labelled tree : " << endl;
   copy_t.display(cout, copy_t.root());
   cout << "... whose structure is then " << endl;
   pt_utree->display(cout, pt_utree->root());
   cout << endl;
   cout << "(" << pt_utree->get_size() << " vertices and depth "
        << pt_utree->get_depth() << ")" << endl;


   // subtree extraction
   cout << "Extracting subtree from this tree rooted at node " << 1 << endl;
   Unlabelled_typed_edge_tree *pt_utreecp= select_subtree(*pt_utree, 1, true);
   pt_utreecp->display(cout, 0);
   delete pt_utreecp;

   // random generation of a tree structure

   cout << "Let's draw a random sample from following distribution " << endl;
   cout << distrib;

   for( int i= 0; i<nv; i++ )
       cout << distrib.simulation() << " ";

   cout << endl;
   cout << "Now a random tree (unlimited depth, " << 2*nv << " vertices max)" << endl;

   // constructor of a Random_unlabelled_tree_generator
   Random_unlabelled_tree_generator<tree_int> generator(distrib,2*nv,2*nv);

   // allocation of an Unlabelled_typed_edge_tree from a Random_unlabelled_tree_generator
   // assignement operator;
   Unlabelled_tree *rd_utree= generator.run();
   Unlabelled_typed_edge_tree *rd_uttree
                                 = new Unlabelled_typed_edge_tree(*rd_utree);

   assert(rd_uttree);

   cout << endl;
   rd_uttree->display(cout, rd_utree->root());

   cout << "... actually : " << rd_uttree->get_size() << " vertices and depth "
        << rd_uttree->get_depth() << endl;

   max_depth= 3;
   cout << "Then a random tree with maximal depth " << max_depth <<
        " and less than " << 3*nv << " vertices" << endl;

   // simulation
   rd_uttree->simulation(var_distrib,3*nv,max_depth);
   cout << "... using following distribution " << endl;
   cout << var_distrib;

   rd_uttree->display(cout, rd_uttree->root());

   // edge types
   cout << "Is (0, 1) an edge ? :" << rd_uttree->is_edge(0, 1) << endl;
   if (rd_uttree->is_edge(0, 1))
      cout << "Type of this edge :" << rd_uttree->edge_type(0, 1) << endl;

   cout << "Is (0, "<< sup_bound+1 <<") an edge ?:"
        << rd_uttree->is_edge(0, sup_bound+1 ) << endl;
   if (rd_uttree->is_edge(0, sup_bound+1))
      cout << "Type of this edge :" << rd_uttree->edge_type(0, sup_bound+1 ) << endl;

   cout << "... actually : " << rd_uttree->get_size() << " vertices and depth "
        << rd_uttree->get_depth() << endl;

   // 2 other simulations
   rd_uttree->simulation(var_distrib,3*nv,max_depth);
   cout << "Another random realization of this process " << endl;
   cout << "(" << rd_uttree->get_size() << " vertices and depth "
        << rd_uttree->get_depth() << ")" << endl;

   rd_uttree->display(cout, rd_uttree->root());

   rd_uttree->simulation(var_distrib,3*nv,max_depth);
   cout << "And another one " << endl;
   cout << "(" << rd_uttree->get_size() << " vertices and depth "
        << rd_uttree->get_depth() << ")" << endl;

   rd_uttree->display(cout, rd_uttree->root());


   // copy constructor of a Unlabelled_tree
   Unlabelled_typed_edge_tree copy_rdtree;
   copy_rdtree= *rd_uttree;
   cout << "A copy from the tree above " << endl;
   copy_rdtree.display(cout, copy_rdtree.root());

   // default Unlabelled_typed_edge_tree
   rd_uttree= new Unlabelled_typed_edge_tree;
   Unlabelled_typed_edge_tree n_utree;

   cout << "Default tree structure is " << endl;
   rd_uttree->display(cout, rd_uttree->root());
   cout << "(" << rd_uttree->get_size() << " vertices and depth "
        << rd_uttree->get_depth() << ")" << endl;

   cout << "or equivalently " << endl;
   n_utree.display(cout, n_utree.root());

   delete rd_uttree;
   delete rd_utree;
   tree_int *tlabellized, *tlabellized2= new tree_int;

   for(i= 0; i < 3; i++)
   {
      // constructor of Unlabelled_typed_edge_tree with no max_depth argument
      rd_uttree= new Unlabelled_typed_edge_tree(distrib, 3*nv);

      cout << "A random structure with less than " << 3*nv << " vertices using no _max_depth constructor " << endl;
      rd_uttree->display(cout, rd_uttree->root());
      cout << "(" << rd_uttree->get_size() << " vertices and depth "
           << rd_uttree->get_depth() << ")" << endl;

      tlabellized= new tree_int(*rd_uttree);
      cout << "Corresponding labellized tree (using constructor): " << endl;
      tlabellized->display(cout, tlabellized->root());
      cout << "(" << tlabellized->get_size() << " vertices and depth "
           << tlabellized->get_depth() << ")" << endl;

      delete tlabellized;

      tlabellized2->set_structure(*rd_uttree, 2);
      cout << "Corresponding labellized tree (using set_structure)" << endl;
      cout << "with default value 2 " << endl;
      tlabellized2->display(cout, tlabellized2->root());
      cout << "(" << tlabellized2->get_size() << " vertices and depth "
           << tlabellized2->get_depth() << ")" << endl;
      delete rd_uttree;
      cout << endl;

   }

   rd_uttree= new Unlabelled_typed_edge_tree(distrib, 3*nv, 3*nv);
   cout << "A random structure with less than " << 3*nv << " vertices using no _max_depth constructor " << endl;
   rd_uttree->display(cout, rd_uttree->root());
   tlabellized= new tree_int(*rd_uttree);

   va= visitor->get_breadthorder(*tlabellized);
   for(i= 0; i < va.size(); i++)
      tlabellized->put(va[i], i+1);

   cout << "Vertices ordered by breadth first traversal: " << endl;
   tlabellized->display(cout, tlabellized->root());

   va= visitor->get_leavesfirstorder(*tlabellized, depth);
   for(i= 0; i < va.size(); i++)
      tlabellized->put(va[i], i+1);

   cout << endl << "Leaves-last tree traversal: " << endl;
   for(i= 0; i < va.size(); i++)
      cout << va[i] << endl;
   cout << endl;

   cout << "Distance of each node to nearest leaf : " << endl;
   for(i= 0; i < va.size(); i++)
      cout << depth[i] << endl;
   cout << endl;

   cout << "Vertices ordered by leaves-first traversal: " << endl;
   tlabellized->display(cout, tlabellized->root());

   cout << "Test of the left (bit) shift operator : " << endl;
   cout << *tlabellized;

   delete rd_uttree;
   delete tlabellized;
   return 0;

}
