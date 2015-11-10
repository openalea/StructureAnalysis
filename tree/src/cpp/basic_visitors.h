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

#include "tree_traits.h"
#include <string>
#include <assert.h>
#include <vector>
#include <deque>
#include <iostream>
#include <cstdlib>


template <class Tree>
struct test_visitor
{
  typedef typename tree_traits<Tree>::vertex_descriptor vid;
  typedef std::vector<vid> vertex_array;

  vertex_array _preorder;
  vertex_array _inorder;
  vertex_array _postorder;

  typedef typename vertex_array::iterator iterator;

  void preorder(vid v, Tree& )
    {
    _preorder.push_back(v);
    }
  void inorder(vid v, Tree& )
    {
    _inorder.push_back(v);
    }
  void postorder(vid v, Tree& )
    {
    _postorder.push_back(v);
    }

  void print( const Tree& t,  std::string order, vertex_array& vertices )
    {
    std::cout << order << ": "<<vertices.size()<<std::endl;

    iterator it= vertices.begin(), end= vertices.end();
    while( it != end )
      std::cout << t.get(*it++) << std::endl;
    std::cout<<std::endl;
    }

  void print_preorder(const Tree& t) { print(t, "preorder", _preorder); }
  void print_inorder(const Tree& t) { print(t, "inorder", _inorder); }
  void print_postorder(const Tree& t) { print(t, "postorder", _postorder); }
};

// mettre en template
// std::string tab= "|-";

template <class Tree>
void display_tree( std::ostream& os,
                   typename tree_traits<Tree>::vertex_descriptor v,
                   Tree& t, std::string tabulation= "")
{
  // mettre en template
  std::string tab= "|-";

  // action
  os<< t.get(v) <<std::endl;
  typename tree_traits<Tree>::children_iterator i, end;
  tie(i, end) = t.children(v);
  if (i != end)
    {
    int s= tabulation.size();
    if( s > 0 )
      {
      tabulation.erase(tabulation.end()-1,tabulation.end());
      tabulation+= " ";
      }
    tabulation+= tab;
    os<<tabulation;

    display_tree(os, *i++, t, tabulation );

    while (i != end)
      {
      os << tabulation;
      display_tree(os, *i++, t, tabulation);
      }
    }
  else
    {
    if( tabulation.size() >= 2 )
      {
      std::string::iterator j=tabulation.end(),i=j-2;
      tabulation.erase(i,j);
      }
    }
}

template<class Tree>
class random_tree
{
public:

  typedef typename tree_traits<Tree>::vertex_descriptor vid;
  typedef typename tree_traits<Tree>::children_iterator ch_it;

  int _max_children;
  int _capacity;
  std::deque<vid> _leaves;
  Tree* _tree;

  random_tree( int max_nb_chilren= 3, int min_max_vertices= 100 ) :
    _max_children(max_nb_chilren),
    _capacity(min_max_vertices),
    _leaves(),
    _tree(0)
    {}

  ~random_tree( )
    {
    _max_children= 0;
    _capacity= 0;
    _tree= 0;
    }

  int run ()
    {
    _tree= new Tree();

    vid root= _tree->add_vertex();
    _tree->root()= root;
    _leaves.push_front(root);
    while( _capacity > 0 )
      random_children(_leaves.back());

    return 0;
    }


  void random_children(vid parent)
    {
    _leaves.pop_back();
    int n= rand() % (_max_children + 1);
    for( int i= 0; i < n; i++ )
      {
      vid v= _tree->add_vertex();
      _tree->add_edge(parent, v);
      _leaves.push_front(v);
      }
    _capacity-= n;
    }

  Tree* tree() { return _tree; }
};
