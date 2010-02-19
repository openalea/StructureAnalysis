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
 *       $Id: int_trees.h 3186 2007-05-25 15:10:30Z dufourko $
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

#ifndef INT_TREES_H
#define INT_TREES_H

namespace Stat_trees
{

/****************************************************************
 *
 *  Purpose :
 *     provides the class "Int_trees" used for representing a
 *     set of oriented trees assumed to be the realizations of an
 *     integral multidimensional tree process
 *     (purpose : iid simulation)
 */

class Int_fl_container;

/****************************************************************
 *
 *  Type definitions
 */

/****************************************************************
 *
 *  Class definitions :
 */

class Int_trees : public Typed_edge_trees<Int_fl_container>
{ // set of tree-structured data assumed to be the realizations
  // of a common inb_integral-dimensional tree process

  // friend classes

// protected :

   // void copy(const Int_trees &otrees);
   // copy operator

public :

   typedef Typed_edge_int_fl_tree<Int_fl_container> tree_type;
   typedef tree_type::value value;

   typedef tree_type::key key;
   typedef tree_type::vertices_size_type vertices_size_type;
   typedef tree_type::children_iterator children_iterator;
   typedef tree_type::vertex_iterator vertex_iterator;

   typedef TreeCharacteristics* TreeCharacteristics_array;
   typedef tree_type** pt_tree_type_array;

   typedef Typed_edge_trees<Int_fl_container> int_trees;
   typedef Typed_edge_trees<Int_fl_container>** pt_int_trees_array;

   Int_trees(int inb_integral = 1, int inb_trees = 0);
   // Int_trees(const Typed_edge_trees<Int_fl_container>& otrees);
   Int_trees(int inb_trees,
                      int* itype,
                      pt_tree_type_array otrees);
   Int_trees(int inb_integral,
             const FrequencyDistribution& ihsize,
             const FrequencyDistribution& ihnb_children,
             bool no_child_flag= false,
             bool init_flag= true);

   virtual ~Int_trees();

   void iid_simulation(const Distribution& distrib);

};
};
#endif
