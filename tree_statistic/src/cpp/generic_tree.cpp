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
 *       $Id: generic_tree.cpp 3193 2007-05-29 10:03:19Z dufourko $
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
#include "generic_tree.h"

using namespace Stat_trees;

/*--------------------------------------------------------------*
 *
 *  Default constructor of Unlabelled_tree class
 *
 *--------------------------------------------------------------*/

Unlabelled_tree::Unlabelled_tree()
 : Generic_tree<char>(0, 1, C_DEFAULT_CHAR)
{ Generic_tree<char>::add_vertex(C_DEFAULT_CHAR); }

/*--------------------------------------------------------------*
 *
 *  Constructor of Unlabelled_tree class using the root id and
 *  the number of vertices
 *
 *--------------------------------------------------------------*/

Unlabelled_tree::Unlabelled_tree(key root,
                                 int n)
 : Generic_tree<char>(root, n, Unlabelled_tree::value())
{}

/*--------------------------------------------------------------*
 *
 *  Copy constructor of Unlabelled_tree class
 *
 *--------------------------------------------------------------*/
Unlabelled_tree::Unlabelled_tree(const Unlabelled_tree& utree)
{ copy(utree); }

/*--------------------------------------------------------------*
 *
 *  Destructor of Unlabelled_tree class
 *
 *--------------------------------------------------------------*/

Unlabelled_tree::~Unlabelled_tree()
{}

/*--------------------------------------------------------------*
 *
 *  Indicate whether the given couple defines an existing edge
 *
 *--------------------------------------------------------------*/


bool Unlabelled_tree::is_edge(key parent, key child)
{ return Generic_tree<char>::is_edge(parent, child); }


/*****************************************************************
*
*  Deletion operator of Unlabelled_tree class
*
**/

void Unlabelled_tree::remove()
{}

/*****************************************************************
 *
 *  Simulation of Unlabelled_tree using a given distribution
 *
 **/

void Unlabelled_tree::simulation(const Distribution& distrib,
                                 int max_size,
                                 int max_depth)

{
   Random_unlabelled_tree_generator< Generic_tree<char> > *tmpg= NULL;
   Unlabelled_tree *tmpt= NULL;
   tmpg=
      new Random_unlabelled_tree_generator< Generic_tree<char> >(distrib, max_size, max_depth);

   tmpt= tmpg->run();
   copy(*tmpt);
   delete tmpg;
   delete tmpt;
   tmpg= NULL;
   tmpt= NULL;
}

/*--------------------------------------------------------------*
 *
 * Copy operator of Unlabelled_tree class
 *
 *--------------------------------------------------------------*/


void Unlabelled_tree::copy(const Unlabelled_tree& utree)
{
  _root= utree._root;
  _children= utree._children;
  _parent= utree._parent;
  _vertices= utree._vertices;

  _size= utree._size;
  _capacity= utree._capacity;

}
