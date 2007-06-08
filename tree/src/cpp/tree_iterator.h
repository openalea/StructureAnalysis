/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR AMAP
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

#ifndef MG_TREE_ITERATOR_H
#define MG_TREE_ITERATOR_H

// Definition de tous les iterateurs sur un MTG
// on peut imaginer:
//
// iterer sur les peres topologiques a partir d'un vertex donne
// iterer sur les peres complexes a partir d'un vertex donne
// iterer sur le contexte topologique a partir d'un vertex donne
// parcourir le mtg a un niveau donne en profondeur d'abord a partir d'un vertex de ce niveau
// idem mais parcours sur les composants d'un vertex
// idem mais parcours sur les contextes topologiques
// parcours en largeur d'abord limite par un parametre (rang, ordre, ...)
//
// pour tous ces iterateurs
//
// - donner la possibilite de parametrer la recherche par des bornes
// - donner la possibilite d'executer une fonction de selection booleenne
// qui permet a l'iterateur de filtrer les elements sur lesquels il itere

#include <iterator>
/*
struct _Rb_tree_base_iterator
{
  typedef _Rb_tree_node_base::_Base_ptr _Base_ptr;
  typedef bidirectional_iterator_tag iterator_category;
  typedef ptrdiff_t difference_type;
  _Base_ptr _M_node;

  void _M_increment()
  {
    if (_M_node->_M_right != 0) {
      _M_node = _M_node->_M_right;
      while (_M_node->_M_left != 0)
        _M_node = _M_node->_M_left;
    }
    else {
      _Base_ptr __y = _M_node->_M_parent;
      while (_M_node == __y->_M_right) {
        _M_node = __y;
        __y = __y->_M_parent;
      }
      if (_M_node->_M_right != __y)
        _M_node = __y;
    }
  }

  void _M_decrement()
  {
    if (_M_node->_M_color == _S_rb_tree_red &&
        _M_node->_M_parent->_M_parent == _M_node)
      _M_node = _M_node->_M_right;
    else if (_M_node->_M_left != 0) {
      _Base_ptr __y = _M_node->_M_left;
      while (__y->_M_right != 0)
        __y = __y->_M_right;
      _M_node = __y;
    }
    else {
      _Base_ptr __y = _M_node->_M_parent;
      while (_M_node == __y->_M_left) {
        _M_node = __y;
        __y = __y->_M_parent;
      }
      _M_node = __y;
    }
  }
};

template < class _Value, class _Ref, class _Ptr >
struct tree_iterator_base
  : public iterator< bidirectional_iterator_tag, _Value,
                     std::ptrdiff_t, Ref, Ptr >
{
  typedef _Value value_type;
  typedef _Ref reference;
  typedef _Ptr pointer;
  typedef tree_iterator_base<_Value, _Value&, _Value*>
    iterator;
  typedef tree_iterator_base<_Value, const _Value&, const _Value*>
    const_iterator;
  typedef tree_iterator_base<_Value, _Ref, _Ptr> _Self;


}
*/

  using namespace std;
  template < typename T >
  struct counting_iterator
    {
    typedef bidirectional_iterator_tag iterator_category;
    typedef T value_type;
    typedef T reference;
    typedef const T* pointer;
    typedef counting_iterator self;

    T _pos;
    counting_iterator( ) : _pos(T(0)) {}
    counting_iterator( T i ): _pos(i) {}

    reference operator*() const { return _pos; }
    pointer operator->() const { return &_pos; }

    self operator++()
      {
      _pos++;
      return *this;
      }
    self operator++(int)
      {
      self tmp= *this;
      _pos++;
      return tmp;
      }
    self operator--()
      {
      _pos--;
      return *this;
      }
    self operator--(int)
      {
      self tmp= *this;
      _pos--;
      return tmp;
      }
    self operator+ (int i) const
      { return self (_pos+i); }
    self operator+= (int i)
      { _pos+= i; return *this; }
    self operator- (int i) const
      { return self (_pos-i); }
    self operator-= (int i)
      { _pos-= i; return *this; }
    self operator- ( self it ) const
      { return self(_pos - it._pos); }
    self operator[] (int i)
      { return self(_pos+i); }
    bool operator < ( self it )
      { return _pos < it._pos; }
    bool operator > ( self it )
      { return _pos > it._pos; }
    bool operator >= ( self it )
      { return _pos >= it._pos; }
    bool operator <= ( self it )
      { return _pos <= it._pos; }
    bool operator == ( self it )
      { return _pos == it._pos; }
    bool operator != ( self it )
      { return _pos != it._pos; }
    };

#endif
// MG_TREE_ITERATOR_H
