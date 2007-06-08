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

#include "tree_simple.h"
#include "tree_traits.h"
#include "tie.h"
#include "basic_visitors.h"
#include <string>
#include <assert.h>
#include <vector>

typedef tree_simple<int> tree;

int main(void)
{
  const int n = 5;

  tree t(0), dt(0);
  typedef tree::vertex_descriptor vertex;
  vertex* v= new vertex[n];
  for( int i= 0; i < n; i++ )
    v[i]= t.add_vertex(i+1);

  // edges
  t.add_edge(v[0], v[1]);
  t.add_edge(v[0], v[2]);
  t.add_edge(v[1], v[3]);
  t.add_edge(v[1], v[4]);

  assert(t.is_root(v[0]));


  vertex root= t.root();

  test_visitor<tree> visitor;
  traverse_tree(t.root(),t,visitor);

  visitor.print_preorder(t);
  visitor.print_inorder(t);
  visitor.print_postorder(t);

  random_tree<tree> generator(3,200);
  generator.run();
  tree* rd_tree= generator.tree();

  assert(rd_tree);

  tree::vertex_iterator it, end;
  tie(it,end)= rd_tree->vertices();
  while( it != end )
    rd_tree->put(*it++,*it+1);

  display_tree(cout, rd_tree->root(), *rd_tree);
  delete rd_tree;

  cout << endl;
  display_tree(cout, t.root(), t);
  cout << endl;

  dt.add_vertex(0);
  display_tree(cout, dt.root(), dt);

  return 0;
}



