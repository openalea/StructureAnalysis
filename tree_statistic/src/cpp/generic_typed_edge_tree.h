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
 *       $Id: generic_typed_edge_tree.h 3193 2007-05-29 10:03:19Z dufourko $
 *
 *       Forum for OpenAlea developers: Openalea-devlp@lists.gforge.inria.f
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

#ifndef GENERIC_TYPED_EDGE_TREE_H
#define GENERIC_TYPED_EDGE_TREE_H

#include "stat_tool/stat_tools.h"
#include "generic_tree.h"

namespace Stat_trees
{
/*! \file generic_typed_edge_tree.h
    \brief Purpose :
           provides the template class "Typed_edge_tree", used for
           the representation of oriented trees with typed edges
*/

/****************************************************************
 *
 *  Class definitions :
 */
class Unlabelled_typed_edge_tree;

/**
   \class Typed_edge_tree
   \brief specialization of the Generic_tree class for
          handling trees with (binary-)typed edges
*/
template <typename Type>
class Typed_edge_tree : public Generic_tree<Type>
{
   template <typename Dummytype> friend
      std::ostream& operator<<(std::ostream& os, Typed_edge_tree<Dummytype>& tree);
   // should be const Typed_edge_tree<Dummytype>& tree
   // friend class Unlabelled_typed_edge_tree;
   template <typename Dummytype> friend class Typed_edge_tree;
   friend class Unlabelled_typed_edge_tree;

public :


   typedef typename Generic_tree<Type>::tree_type tree_type;
   typedef typename tree_type::value value;

   typedef typename tree_type::key key;
   typedef typename tree_type::vertices_size_type vertices_size_type;
   typedef typename tree_type::children_iterator children_iterator;
   typedef typename tree_type::vertex_iterator vertex_iterator;

   /// false is intended to represent branching +
   /// and true to represent succession <
   typedef std::vector<bool> type_list;


protected :

   /// liste of edge types
   type_list edge_types;

   void display_skip(ostream& os,
                     key v,
                     std::string tabulation,
                     bool& lineskip);

   void copy(const Typed_edge_tree<Type>& tree);

public :

   Typed_edge_tree(key root= 0,
                   int n= 1,
                   const value& default_value= value());
   Typed_edge_tree(const Typed_edge_tree<Type>& tree);
   Typed_edge_tree(const Generic_tree<Type>& tree);
   Typed_edge_tree(Unlabelled_typed_edge_tree& utree,
                   const value& default_value= value());
   ~Typed_edge_tree();
   Typed_edge_tree<Type>& operator=(const Typed_edge_tree<Type>& tree);

   /** return the structure of \e self */
   Unlabelled_typed_edge_tree* get_structure(); // const
   /** creation of a tree with given structure and default attribute values */
   void set_structure(Unlabelled_typed_edge_tree& utree,
                      const value& default_value= value());
   /** add a vertex and return its identifier */
   key add_vertex(const value& data= value());
   /** create a typed edge or set the type of an existing edge */
   bool add_edge(key parent, key child, bool type= false);
   /** return the type of a given edge */
   bool edge_type(key parent, key child);
   /** set the type of a given edge */
   bool set_edge_type(key parent, key child, bool type= false);
   /** return the branching order (depth) of the tree */
   unsigned int get_branching_order();
   /** return the branching order (depth) of a given vertex */
   unsigned int get_branching_order(key v);

};

template <class Tree>
Tree* select_typed_edge_subtree(Tree& t,
                                typename tree_traits<Tree>::vertex_descriptor v,
                                bool keep_flag= true);

typedef Typed_edge_tree<char> Generic_typed_edge_tree_char;

class Unlabelled_typed_edge_tree : public Generic_typed_edge_tree_char,
                                   public Unlabelled_tree
{  // a tree structure, i.e. a tree with no label

   // friend classes
   // friend  class Typed_edge_tree<char>;

public :


   typedef Generic_typed_edge_tree_char::tree_type tree_type;
   // typedef Unlabelled_tree::tree_type tree_type;
   typedef Unlabelled_tree::value value;

   typedef tree_type::key key;
   typedef tree_type::vertices_size_type vertices_size_type;
   typedef tree_type::children_iterator children_iterator;
   typedef tree_type::vertex_iterator vertex_iterator;

   /// false is intended to represent branching +
   /// and true to represent succession <
   /// typedef std::vector<bool> type_list;

protected:

   void copy(const Unlabelled_typed_edge_tree &utree);
   void copy(const Unlabelled_tree &utree);
   void remove();

public:

   Unlabelled_typed_edge_tree();
   Unlabelled_typed_edge_tree(key root,
                              int n,
                              const value& default_value= value());
   Unlabelled_typed_edge_tree(const Unlabelled_typed_edge_tree& utree);
   /** Return the structure of a tree */
   template<typename Type> Unlabelled_typed_edge_tree(Generic_tree<Type>& gtree);
   /** Return the structure of a typed-edge tree */
   template<typename Type> Unlabelled_typed_edge_tree(Typed_edge_tree<Type>& gtree);
   Unlabelled_typed_edge_tree(const Distribution& distrib,
                              int max_size= I_DEFAULT_TREE_SIZE,
                              int max_depth= I_DEFAULT_TREE_DEPTH);
   ~Unlabelled_typed_edge_tree();
   Unlabelled_typed_edge_tree& operator=(const Unlabelled_typed_edge_tree &tree);

   /** return root vertex*/
   key root() const;
   /** return whether argument is the root vertex*/
   bool is_root(key vid);
   /** get property */
   const value& get(key vid);
   /** put property  */
   void put(key vid, const value& v= C_DEFAULT_CHAR);
   /** add a vertex and return its identifier */
   key add_vertex();
   key add_vertex(const value& dummy_value);
   /** return parent vertex */
   key parent(key p) const;
   /** return children vertices */
   pair_ch_it children(key v);
   /** return vertices */
   pair_it vertices();
   /** return the tree size */
   unsigned int get_size() const;
   /** return tree depth */
   unsigned int get_depth();
   /** return depth of a vertex */
   unsigned int get_depth(key v);
   /** return true iif the given couple defines an edge */
   bool is_edge(key parent, key child);
   /** create a typed edge or set the type of an existing edge */
   bool add_edge(key parent, key child, bool type= false);

   /** print tree */
   void display(ostream& os, key v);
   void simulation(const Distribution& distrib,
                   int max_size= I_DEFAULT_TREE_SIZE,
                   int max_depth= I_DEFAULT_TREE_DEPTH);
   void simulation(const FrequencyDistribution& hist,
                   int max_size= I_DEFAULT_TREE_SIZE,
                   int max_depth= I_DEFAULT_TREE_DEPTH);

};

#include "generic_typed_edge_tree.hpp"
};

#endif
