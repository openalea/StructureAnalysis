/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): J.-B. Durand (jean-baptiste.durand@imag.fr)
 *
 *       $Source$
 *       $Id: generic_tree.h 3193 2007-05-29 10:03:19Z dufourko $
 *
 *       Forum for V-Plants developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS for A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef GENERIC_TREE_H
#define GENERIC_TREE_H

#include "stat_tool/stat_tools.h"
#include <deque>

using Tree_tie::tie;

namespace Stat_trees
{

/****************************************************************
 *
 *  Purpose :
 *     provides the template class "generic_tree", used for
 *     representing oriented trees with different label types
 */

class Unlabelled_tree;

/*************************************************************
 *
 *  Constants :
 */


const int I_DEFAULT_TREE_SIZE= 50;             // default tree size
const int I_DEFAULT_TREE_DEPTH= 5;             // default tree depth
const int I_DEFAULT_INF_BOUND= 0;              // default minumum number of children
const int I_DEFAULT_SUP_BOUND= 2;              // default maximum number of children
const int I_DEFAULT_MAX_CHILDREN= I_DEFAULT_SUP_BOUND;
const double D_DEFAULT_PROBABILITY= 0.5;       // parameter value for default distribution
const int I_DEFAULT_IDENT= BINOMIAL;           // default type of distribution
const double D_DEFAULT_PARAMETER= D_DEFAULT;   // cf. Distribution()
const char C_DEFAULT_CHAR= '*';                // default label for char trees
const DiscreteParametric P_DEFAULT_DISTRIBUTION(I_DEFAULT_SUP_BOUND-I_DEFAULT_INF_BOUND+1,
                                                I_DEFAULT_IDENT,
                                                I_DEFAULT_INF_BOUND,
                                                I_DEFAULT_SUP_BOUND,
                                                D_DEFAULT_PARAMETER,
                                                D_DEFAULT_PROBABILITY);
                                               // default distribution

/****************************************************************
 *
 *  Class definitions :
 */

template <class Tree>
struct generic_visitor
{
   typedef typename tree_traits<Tree>::vertex_descriptor vid;
   typedef std::vector<vid> vertex_array;
   typedef std::deque<vid> vertex_deque;

   vertex_array _preorder;
   vertex_array _inorder;
   vertex_array _postorder;

   typedef typename vertex_array::iterator iterator;

   void preorder(vid v, Tree& )
     { _preorder.push_back(v); }

   void inorder(vid v, Tree& )
     { _inorder.push_back(v); }

   void postorder(vid v, Tree& )
     { _postorder.push_back(v); }

   void print( const Tree& t,  std::string order, vertex_array& vertices )
   {
      cout << order << ": "<<vertices.size()<<endl;

      iterator it= vertices.begin(), end= vertices.end();
      while( it != end )
        cout << t.get(*it++) << endl;
      cout<<endl;
   }

   void print_preorder(const Tree& t) { print(t, "preorder", _preorder); }
   void print_inorder(const Tree& t) { print(t, "inorder", _inorder); }
   void print_postorder(const Tree& t) { print(t, "postorder", _postorder); }

   vertex_array get_preorder(const Tree& t) { return _preorder; }
   vertex_array get_inorder(const Tree& t) { return _inorder; }
   vertex_array get_postorder(const Tree& t) { return _postorder; }
   vertex_array get_breadthorder(Tree& t, vid vroot) // should be const Tree& t
   {   // breadth-first tree traversing

       typedef typename tree_traits<Tree>::children_iterator children_iterator;
       int curr_width, nb_children;
       std::vector<vid> res(0);
       std::deque<vid> nodes(0);
       vid v;
       children_iterator it, end;

       curr_width= 1; // number of remaining children at current width
       nb_children= 0; // number of children at next width
       nodes.push_front(vroot);

       while(!nodes.empty())
       // breadth-first tree traversing
       {
          if (!curr_width)
            {
              curr_width= nb_children;
              nb_children= 0;
            }

           v= nodes.back();
           nodes.pop_back();
           res.push_back(v);
           curr_width-= 1;
           // current vertex need not be processed any more

           Tree_tie::tie(it, end)= t.children(v);

           while (it < end)
           {
              nb_children++;
              nodes.push_front(*it++);
           }
       }
       return res;
   }
   vertex_array get_breadthorder(Tree& t) // should be const Tree& t
   {   // breadth-first tree traversing from root
      std::vector<vid> res= get_breadthorder(t, t.root());

      assert(res.size() == (unsigned int)t.size());
      return res;
   }

   void find_leaves(Tree& t,
                    vid current_node,
                    std::vector<vid>& res)
   {
      typedef typename tree_traits<Tree>::children_iterator children_iterator;
      children_iterator it, end;

      Tree_tie::tie(it, end)= t.children(current_node);
      if (it == end)
         res.push_back(current_node);
      else
         while (it < end)
            find_leaves(t, *it++, res);
   }

   vertex_array get_leavesfirstorder(Tree& t, std::vector<int>& depth) // should be const Tree& t
   {   // tree-traversing where the leaf nodes come first, then their parents, etc.
       // the distance of each node to nearest leaf is given by the table "depth"

       typedef typename tree_traits<Tree>::children_iterator children_iterator;
       typedef typename tree_traits<Tree>::children_iterator children_iterator;

       const unsigned int size= t.get_size();
       unsigned int i;
       std::deque<vid> inner_nodes;
       std::vector<vid> res;
       unsigned int curr_width, next_width, curr_depth;
       vid p, current_child;
       children_iterator it, end;
       bool* occurred; // bool[i] == true <=> i exists in "res"


       find_leaves(t, t.root(), res);
       occurred= new bool[size];
       depth.resize(0);

       for(i= 0; i < size; i++)
          occurred[i]= false;

       curr_width= res.size();
       // number of remaining nodes before the depth
       // increases (starting from the leaf nodes)
       next_width= 0;
       // number of nodes having constant depth == curr_depth+1
       curr_depth= 0;

       for(i= 0; i < curr_width; i++)
       {
          inner_nodes.push_back(res[i]);
          occurred[res[i]]= true;
          depth.push_back(0); // distance from a leaf node to itself
       }

       while (res.size() < size)
       {
          current_child= inner_nodes.back();
          if (!curr_width)
          {
             curr_width= next_width;
             next_width= 0;
             curr_depth++;
          }

          if (current_child != t.root())
          {
             p= t.parent(current_child);
             inner_nodes.push_front(p);
             next_width++;
             if (!occurred[p])
             {
                occurred[p]= true;
                depth.push_back(curr_depth+1);
                res.push_back(p);
             }
          }
          inner_nodes.pop_back();
          curr_width--;
       }

       assert(res.size() == t.get_size());
       assert(depth.size() == t.get_size());
       delete [] occurred;
       occurred= NULL;

       return res;
   }

   /** Return the array of all ancestors of a given vertex (including this vertex),
       beginning by root vertex */
   vertex_deque get_vertex_ancestors(Tree& t,
                                     typename tree_traits<Tree>::vertex_descriptor v) // should be const Tree& t
   {   typedef typename tree_traits<Tree>::vertex_descriptor vertex_descriptor;

       vertex_descriptor current= v;
       vertex_deque res;

       res.push_front(current);
       while(current != t.root())
       {
          current= t.parent(current);
          res.push_front(current);
       }
       return res;
   }
};

template <class Tree>
Tree* select_subtree(Tree& t,
                     typename tree_traits<Tree>::vertex_descriptor v,
                     bool keep_flag= true);

template <typename Type>
class Generic_tree : public tree_simple<Type>
{  // class for the representation of tree-structured data with
   // labelled vertices of a generic type

   template <typename Dummytype> friend
      std::ostream& operator<<(std::ostream& os, Generic_tree<Dummytype>& tree);
   // should be const Generic_tree<Dummytype>& tree
   // friend class Unlabelled_tree;

public :

   typedef tree_simple<Type> tree_type;
   typedef typename tree_type::value value;

   typedef typename tree_type::key key;
   typedef typename tree_type::vertices_size_type vertices_size_type;
   typedef typename tree_type::children_iterator children_iterator;
   typedef typename tree_type::vertex_iterator vertex_iterator;

protected :

   void display_skip(ostream& os,
                     key v,
                     std::string tabulation,
                     bool& lineskip);

   void copy(const Generic_tree<Type>& tree);

public :

   Generic_tree(key root= 0,
                int n= 1,
                const value& default_value= value());
   Generic_tree(const Generic_tree<Type>& tree);
   Generic_tree(Unlabelled_tree& utree, const value& default_value= value());
   ~Generic_tree();
   Generic_tree<Type>& operator=(const Generic_tree<Type>& tree);

   // deletion of the tree labels
   Unlabelled_tree* get_structure(); // const
   void set_structure(Unlabelled_tree& utree,
                      const value& default_value= value());

   unsigned int get_size() const;
   unsigned int get_depth();
   unsigned int get_depth(key v);
   unsigned int get_nb_children(key v);

   /** return true iif the given couple defines an edge */
   bool is_edge(key parent, key child);
   void display(ostream& os, key v); // should be const
};

typedef Generic_tree<char> Generic_tree_char;

class Unlabelled_tree : public Generic_tree_char
{  // a tree structure, i.e. a tree with no label

   // friend classes
   // friend class Generic_tree<char>;

protected:

   void copy(const Unlabelled_tree &utree);
   void remove();

public:

   Unlabelled_tree();
   Unlabelled_tree(key root,
                   int n);
   Unlabelled_tree(const Unlabelled_tree& utree);
   /** Return the structure of a tree */
   template<typename Type> Unlabelled_tree(Generic_tree<Type>& gtree);
   Unlabelled_tree(const Distribution& distrib,
                   int max_size= I_DEFAULT_TREE_SIZE,
                   int max_depth= I_DEFAULT_TREE_DEPTH);
   ~Unlabelled_tree();
   Unlabelled_tree& operator=(const Unlabelled_tree &tree);

   void simulation(const Distribution& distrib,
                   int max_size= I_DEFAULT_TREE_SIZE,
                   int max_depth= I_DEFAULT_TREE_DEPTH);
   void simulation(const FrequencyDistribution& hist,
                   int max_size= I_DEFAULT_TREE_SIZE,
                   int max_depth= I_DEFAULT_TREE_DEPTH);

   /** return true iif the given couple defines an edge */
   bool is_edge(key parent, key child);

};

template<class Tree>
class Random_unlabelled_tree_generator
{  // random generator of a tree structure

private:

   int _curr_width; // size of the unprocessed forest
   int _nb_children; // number of spawned children

public:

   typedef typename tree_traits<Tree>::vertex_descriptor vid;
   typedef typename tree_traits<Tree>::children_iterator ch_it;

   Distribution* _distribution;
   int _max_size;
   int _max_depth;
   std::deque<vid> _leaves;
   // deque of the next leaves to process

   Random_unlabelled_tree_generator(const Distribution& distrib= P_DEFAULT_DISTRIBUTION,
                                    int max_size= I_DEFAULT_TREE_SIZE,
                                    int max_depth= I_DEFAULT_TREE_DEPTH);

   Random_unlabelled_tree_generator(const Random_unlabelled_tree_generator& rutg);

   ~Random_unlabelled_tree_generator();

   Unlabelled_tree*  run();
   void random_children(vid parent, Unlabelled_tree* const _tree);
};


#include "generic_tree.hpp"
};

#endif
