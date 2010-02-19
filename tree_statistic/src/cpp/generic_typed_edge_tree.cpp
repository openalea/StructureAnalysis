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
 *       $Id: generic_typed_edge_tree.cpp 3193 2007-05-29 10:03:19Z dufourko $
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

#include "tree/tree_simple.h"
#include "tree/tree_traits.h"
#include "tree/basic_visitors.h"
#include "stat_tool/stat_tools.h"
#include "generic_typed_edge_tree.h"
#include <deque>

using namespace Stat_trees;

/*****************************************************************
 *
 *  Default constructor of Unlabelled_tree class
 *
 **/

Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree()
 : Typed_edge_tree<char>(0, 1, C_DEFAULT_CHAR)
 , Unlabelled_tree()
{ Typed_edge_tree<char>::add_vertex(C_DEFAULT_CHAR); }

/*****************************************************************
 *
 *  Constructor of Unlabelled_typed_edge_tree class specifying root and
 *  number of vertices
 *
 **/

Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree(key root,
                                                       int n,
                                                       const value& default_value)
 : Typed_edge_tree<char>(root, n, default_value)
 , Unlabelled_tree(root, n)
{
  Unlabelled_tree::_root= Typed_edge_tree<char>::_root;
  Unlabelled_tree::_children= Typed_edge_tree<char>::_children;
  Unlabelled_tree::_parent= Typed_edge_tree<char>::_parent;
  Unlabelled_tree::_vertices= Typed_edge_tree<char>::_vertices;

  Unlabelled_tree::_size= Typed_edge_tree<char>::_size;
  Unlabelled_tree::_capacity= Typed_edge_tree<char>::_capacity;
}

/*****************************************************************
 *
 *  Copy constructor of Unlabelled_typed_edge_tree class
 *
 **/

Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree(const Unlabelled_typed_edge_tree& utree)
 : Typed_edge_tree<char>(utree)
 , Unlabelled_tree(utree)
{ copy(utree); }

/*****************************************************************
 *
 *  Constructor of Unlabelled_typed_edge_tree class using a distribution
 *
 **/

Unlabelled_typed_edge_tree::Unlabelled_typed_edge_tree(const Distribution& distrib,
                                                       int max_size,
                                                       int max_depth)
{
   Typed_edge_tree<char>(0, 1, C_DEFAULT_CHAR);
   simulation(distrib, max_size, max_depth);
}

/*****************************************************************
 *
 *  Destructor of Unlabelled_typed_edge_tree class
 *
 **/

Unlabelled_typed_edge_tree::~Unlabelled_typed_edge_tree()
{}

/*****************************************************************
 *
 *  Deletion operator of Unlabelled_typed_edge_tree class
 *
 **/

void Unlabelled_typed_edge_tree::remove()
{}


/*--------------------------------------------------------*
 *
 *  Assignment operator of Unlabelled_typed_edge_tree class
 *
 **/


Unlabelled_typed_edge_tree& Unlabelled_typed_edge_tree::operator=(const Unlabelled_typed_edge_tree& utree)
{
  if (&utree != this)
    {
      remove();
      copy(utree);
    }

  return *this;
}
/*--------------------------------------------------------------*
 *
 * Return the root vertex of an Unlabelled_typed_edge_tree
 *
 *--------------------------------------------------------------*/

Unlabelled_typed_edge_tree::key Unlabelled_typed_edge_tree::root() const
{ return Unlabelled_tree::root(); }

/*--------------------------------------------------------------*
 *
 * Return whether the argument is the root vertex of an Unlabelled_typed_edge_tree
 *
 *--------------------------------------------------------------*/

bool Unlabelled_typed_edge_tree::is_root(key vid)
{ return Unlabelled_tree::is_root(vid); }

/*--------------------------------------------------------------*
 *
 * Return property associated with given vertex for Unlabelled_typed_edge_tree
 *
 *--------------------------------------------------------------*/

const Unlabelled_typed_edge_tree::value& Unlabelled_typed_edge_tree::get(key vid)
{ return Unlabelled_tree::get(vid); }


/*--------------------------------------------------------------*
 *
 * Set property associated with given vertex for Unlabelled_typed_edge_tree
 *
 *--------------------------------------------------------------*/

void Unlabelled_typed_edge_tree::put(key vid, const value& v)
{ return Unlabelled_tree::put(vid, v); }

/*--------------------------------------------------------------*
 *
 * Return the size of an Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

unsigned int Unlabelled_typed_edge_tree::get_size() const
{ return Unlabelled_tree::get_size(); }

/*--------------------------------------------------------------*
 *
 * Return the depth of an Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

unsigned int Unlabelled_typed_edge_tree::get_depth()
{ return Unlabelled_tree::get_depth(); }

/*--------------------------------------------------------------*
 *
 * Return the vertex depth for an Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

unsigned int Unlabelled_typed_edge_tree::get_depth(key v)
{ return Unlabelled_tree::get_depth(v); }

/*--------------------------------------------------------------*
 *
 * Return the parent of a vertex for Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

Unlabelled_typed_edge_tree::key Unlabelled_typed_edge_tree::parent(key p) const
{ return Unlabelled_tree::parent(p); }

/*--------------------------------------------------------------*
 *
 * Return the children of a vertex for Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

Unlabelled_typed_edge_tree::pair_ch_it
Unlabelled_typed_edge_tree::children(key v)
{ return Unlabelled_tree::children(v); }

/*--------------------------------------------------------------*
 *
 * Return the vertices of an Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

Unlabelled_typed_edge_tree::pair_it
Unlabelled_typed_edge_tree::vertices()
{ return Unlabelled_tree::vertices(); }

/*--------------------------------------------------------------*
 *
 * Print an Unlabelled_typed_tree
 *
 *--------------------------------------------------------------*/

void Unlabelled_typed_edge_tree::display(ostream& os, key v)
{ return Unlabelled_tree::display(os, v); }

/*--------------------------------------------------------------*
 *
 * Add a vertex to Unlabelled_typed_tree and return its identifier
 *
 *--------------------------------------------------------------*/

Unlabelled_typed_edge_tree::key Unlabelled_typed_edge_tree::add_vertex()
{
   edge_types.push_back(false);
   Typed_edge_tree<char>::add_vertex(C_DEFAULT_CHAR);
   return Unlabelled_tree::add_vertex(C_DEFAULT_CHAR);
}

Unlabelled_typed_edge_tree::key Unlabelled_typed_edge_tree::add_vertex(const value& dummy_value)
{
   edge_types.push_back(false);
   Typed_edge_tree<char>::add_vertex(C_DEFAULT_CHAR);
   return Unlabelled_tree::add_vertex(C_DEFAULT_CHAR);
}

/*--------------------------------------------------------------*
 *
 *  Return true iif the given couple defines an edge
 *
 *--------------------------------------------------------------*/

bool Unlabelled_typed_edge_tree::is_edge(key parent,
                                         key child)
{ return Unlabelled_tree::is_edge(parent, child); }


/*--------------------------------------------------------------*
 *
 *  Add an edge to an Unlabelled_typed_edge_tree
 *
 *--------------------------------------------------------------*/

bool Unlabelled_typed_edge_tree::add_edge(key parent,
                                          key child,
                                          bool type)
{
   bool res= true;
   if (!Unlabelled_tree::is_edge(parent, child))
   {
      Typed_edge_tree<char>::add_edge(parent, child);
      res= Unlabelled_tree::add_edge(parent, child);
   }
   edge_types[child]= type;
   return res;
}


/*****************************************************************
 *
 *  Simulation of Unlabelled_typed_edge_tree using a given distribution
 *
 **/

void Unlabelled_typed_edge_tree::simulation(const Distribution& distrib,
                                            int max_size,
                                            int max_depth)

{
   Random_unlabelled_tree_generator< Typed_edge_tree<char> > *tmpg= NULL;
   Unlabelled_tree *tmpt= NULL;
   tmpg=
      new Random_unlabelled_tree_generator< Typed_edge_tree<char> >(distrib, max_size, max_depth);

   tmpt= tmpg->run();
   copy(*tmpt);
   delete tmpg;
   delete tmpt;
   tmpg= NULL;
   tmpt= NULL;
}

/*****************************************************************
 *
 *  Simulation of Unlabelled_typed_edge_tree using a given frequency distribution
 *
 **/

void Unlabelled_typed_edge_tree::simulation(const FrequencyDistribution& hist,
                                            int max_size,
                                            int max_depth)
{
   Distribution distrib(hist);

   simulation(distrib, max_size, max_depth);
}

/*****************************************************************
 *
 *  Copy operator of Unlabelled_typed_edge_tree class
 *
 **/


void Unlabelled_typed_edge_tree::copy(const Unlabelled_typed_edge_tree& utree)
{
  Typed_edge_tree<char>::_root= utree.Typed_edge_tree<char>::_root;
  Typed_edge_tree<char>::_children= utree.Typed_edge_tree<char>::_children;
  Typed_edge_tree<char>::_parent= utree.Typed_edge_tree<char>::_parent;
  Typed_edge_tree<char>::_vertices= utree.Typed_edge_tree<char>::_vertices;

  Typed_edge_tree<char>::_size= utree.Typed_edge_tree<char>::_size;
  Typed_edge_tree<char>::_capacity= utree.Typed_edge_tree<char>::_capacity;

  Unlabelled_tree::_root= Typed_edge_tree<char>::_root;
  Unlabelled_tree::_children= Typed_edge_tree<char>::_children;
  Unlabelled_tree::_parent= Typed_edge_tree<char>::_parent;
  Unlabelled_tree::_vertices= Typed_edge_tree<char>::_vertices;

  Unlabelled_tree::_size= Typed_edge_tree<char>::_size;
  Unlabelled_tree::_capacity= Typed_edge_tree<char>::_capacity;
}

/*****************************************************************
 *
 *  Copy operator of Unlabelled_typed_edge_tree class
 *
 **/


void Unlabelled_typed_edge_tree::copy(const Unlabelled_tree& utree)
{
  Typed_edge_tree<char>::_root= utree._root;
  Typed_edge_tree<char>::_children= utree._children;
  Typed_edge_tree<char>::_parent= utree._parent;
  Typed_edge_tree<char>::_vertices= utree._vertices;

  Typed_edge_tree<char>::_size= utree.size();
  // Typed_edge_tree<char>::_capacity= utree._capacity;
  Typed_edge_tree<char>::edge_types.resize(utree.size(), false);

  Unlabelled_tree::_root= Typed_edge_tree<char>::_root;
  Unlabelled_tree::_children= Typed_edge_tree<char>::_children;
  Unlabelled_tree::_parent= Typed_edge_tree<char>::_parent;
  Unlabelled_tree::_vertices= Typed_edge_tree<char>::_vertices;

  Unlabelled_tree::_size= Typed_edge_tree<char>::_size;
  // Unlabelled_tree::_capacity= Typed_edge_tree<char>::_capacity;

}
