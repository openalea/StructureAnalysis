/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): C. Pradal (christophe.pradal@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
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
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */

#ifndef MG_TREE_TRAITS_H
#define MG_TREE_TRAITS_H

#include "tie.h"

// type de reference

template < class Tree >
struct tree_traits
{
  typedef typename Tree::vertex_descriptor  vertex_descriptor;
  typedef typename Tree::children_iterator  children_iterator;

  typedef typename Tree::vertex_iterator  vertex_iterator;
  typedef typename Tree::vertices_size_type vertices_size_type;

//  typedef typename T::edge_descriptor  edge_descriptor;
};

template <class Tree, class TreeVisitor>
void traverse_tree( typename tree_traits<Tree>::vertex_descriptor v,
                    Tree& t, TreeVisitor& visitor )
{
  visitor.preorder(v, t);
  typename tree_traits<Tree>::children_iterator i, end;
  tie(i, end) = t.children(v);
  if (i != end)
    {
    traverse_tree(*i++, t, visitor);
    visitor.inorder(v, t);
    while (i != end)
      traverse_tree(*i++, t, visitor);
    }
  else
    visitor.inorder(v, t);
  visitor.postorder(v, t);
}

struct null_tree_visitor
{
  template <typename Vertex, typename Tree> void preorder(Vertex, Tree&) { }
  template <typename Vertex, typename Tree> void inorder(Vertex, Tree&) { }
  template <typename Vertex, typename Tree> void postorder(Vertex, Tree&) { }
};


#endif
// MG_TREE_TRAITS_H
