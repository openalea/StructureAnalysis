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

#ifndef MG_TREE_SIMPLE_H
#define MG_TREE_SIMPLE_H

#include <utility>
#include <vector>
#include "tree_iterator.h"
#include <assert.h>

/* ----------------------------------------------------------------------- */

/**
   \class vertex_simple
   \brief A container of internal vertex property
*/

template < typename T >
struct vertex_simple
{
  typedef T value;
  value _data;

  vertex_simple( const value& data= value() ) : _data(data) {}

  /// get property
  const value& get() const {return _data;}
  /// put property
  void put( const value& data ) { _data= data; }
};

/* ----------------------------------------------------------------------- */

/**
   \class tree_simple
   \brief A non mutable tree container where vertex are just index.
*/

template < typename T, template< typename > class V = vertex_simple >
class tree_simple
{
public:
  typedef V< T > vertex;
  typedef int key;
  typedef typename vertex::value value;

  typedef key vertex_descriptor;
  typedef int vertices_size_type;
  typedef std::vector<key>::iterator children_iterator;
  typedef counting_iterator<key> vertex_iterator;

  typedef std::pair<children_iterator,children_iterator> pair_ch_it;
  typedef std::pair<vertex_iterator,vertex_iterator> pair_it;

public:
  key _root;
  std::vector< std::vector< key > > _children;
  std::vector< key > _parent;
  std::vector< vertex > _vertices;

protected:
  unsigned int _size; // int size;
  int _capacity;

public:
  tree_simple( key root= 0,
               int n= 1,
               const value& default_value= value() ):
    _root(root),
    _children(),
    _parent(),
    _vertices(),
    _size(0)
    {
    _children.reserve(n);
    _parent.reserve(n);
    _vertices.reserve(n);

//    for(int i= 0; i < n; i++)
//      add_vertex();

    }
  tree_simple( const tree_simple& t ):
    _root(t._root),
    _children(t._children),
    _parent(t._parent),
    _vertices(t._vertices),
    _size(t._size)
    {}
  ~tree_simple()
    {}

   vertex_descriptor root() const
    { return _root; }

   vertex_descriptor& root()
    { return _root; }

   vertex_descriptor parent( vertex_descriptor node ) const
    {
    assert( (unsigned int)node < _size );
    return _parent[node];
    }
  pair_ch_it children( vertex_descriptor node )
    {
    assert( (unsigned int)node < _size );
    return pair_ch_it( _children[node].begin(), _children[node].end() );
    }
  int size() const
    {  return _size; }

  pair_it vertices()
    {
    return pair_it( vertex_iterator(0), vertex_iterator(_size) );
    }

  bool is_root( vertex_descriptor node )
    { return ( node == _root ) || ( _parent[node] == -1 ); }
  bool is_leaf( vertex_descriptor node )
    { return _children[node].size() == 0; }

  vertex_descriptor add_vertex( const value& data= value() )
    {
    _size++;
    _vertices.push_back(vertex(data));
    _parent.push_back(vertex_descriptor(-1));
    std::vector<key> no_child;
    _children.push_back(no_child);

    assert(_vertices.size()==_size);
    assert(_parent.size()==_size);
    assert(_children.size()==_size);

    return _size-1;
    }

/*
  void remove_vertex( vertex_descriptor node )
    {
    _size--;
    _vertices.erase( _vertices.begin()+int(node) );
    _parent.erase( _parent.begin()+int(node) );
    _children.erase( _children.begin()+int(node) );
    assert(_vertices.size()==_size);
    assert(_parent.size()==_size);
    assert(_children.size()==_size);
    }
*/
  bool add_edge( vertex_descriptor parent, vertex_descriptor child )
    {
    assert( (unsigned int)parent < _size );
    assert( (unsigned int)child < _size );
    _parent[child]= parent;
    _children[parent].push_back(child);
    return true;
    }

  /// get property
  const value& get( vertex_descriptor node ) const
    { return _vertices[node].get(); }

  /// put property
  void put( vertex_descriptor node, const value& data )
    { _vertices[node].put(data); }

};



#endif
// MG_TREE_SIMPLE_H
