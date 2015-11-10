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

#ifndef MG_TREE_CONCEPTS_H
#define MG_TREE_CONCEPTS_H

template < class T >
class TreeConcept
{
  typedef typename tree_traits<T>::vertex_descriptor vertex_descriptor;
  typedef typename tree_traits<T>::children_iterator children_iterator;
  typedef typename tree_traits<T>::vertices_size_type vertices_size_type;
  typedef typename tree_traits<T>::vertex_iterator vertex_iterator;

  virtual vertex_descriptor root() = 0;
  virtual vertex_descriptor parent(vertex_descriptor v) = 0;
  virtual std::pair<children_iterator,children_iterator>
    children(vertex_descriptor v) = 0;
}

template < class T >
class MutableTreeConcept : public TreeConcept< T >
{
  virtual vertex_descriptor add_vertex( vertex_descriptor parent ) = 0;
  virtual void remove_vertex( vertex_descriptor u ) = 0;
}

template <class T>
class VertexListTreeConcept : public TreeConcept< T >
{
  typedef typename tree_traits<T>::vertex_iterator vertex_iterator;
  typedef typename tree_traits<T>::vertices_size_type vertices_size_type;
  typedef typename tree_traits<T>::vertices_size_type vertices_size_type;

  virtual vertex_iterator vertices() = 0;
  virtual vertices_size_type size() = 0;
};


#endif
// MG_TREE_CONCEPTS_H
