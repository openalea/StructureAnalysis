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
 *       $Id: int_trees.cpp 3193 2007-05-29 10:03:19Z dufourko $
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
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
// #include "sequence_analysis/sequences.h"

#include "int_fl_containers.h"
#include "generic_typed_edge_tree.h"
#include "typed_edge_trees.h"
#include "int_trees.h"
#include <iostream>
#include <cstdlib>

using namespace Stat_trees;

/*****************************************************************
 *
 *  Default constructor of Observed_int_tree class
 *
 **/

Int_trees::Int_trees(int inb_integral, int inb_trees)
 : Typed_edge_trees<Int_fl_container>(inb_integral, 0, inb_trees)
{}


/*****************************************************************
 *
 *  Constructor of Int_trees class using the number of trees,
 *  the type of each variable and multidimensional observed trees
 *
 **/

Int_trees::Int_trees(int inb_trees,
                     int* itype,
                     Typed_edge_int_fl_tree<Int_fl_container >** otrees)
{
   value v;
   Typed_edge_trees<Int_fl_container> *tmp;
   if (otrees != NULL)
   {
        v= (*otrees)[0].get((*otrees)[0].root());
        assert(v.nb_float() == 0);
        tmp= new Typed_edge_trees<Int_fl_container>(inb_trees, itype, otrees);
        Typed_edge_trees<Int_fl_container>::copy(*tmp);
        delete tmp;
        tmp= NULL;
   }
}

/*****************************************************************
 *
 *  Constructor of Int_trees class
 *  using the number of integral variables,
 *  frequency distributions for the size and number of children of the trees,
 *  a flag on the possibility for a node to have no child due to random
 *  and a flag on the (random) tree initialization
 *
 **/

Int_trees::Int_trees(int inb_integral,
                     const FrequencyDistribution& ihsize,
                     const FrequencyDistribution& ihnb_children,
                     bool no_child_flag,
                     bool init_flag)
 : Typed_edge_trees<Int_fl_container>(inb_integral, 0,
                                    ihsize, ihnb_children,
                                    no_child_flag, init_flag)
{}

/*****************************************************************
 *
 *  Destructor of Int_trees class
 *
 **/

Int_trees::~Int_trees()
{}

/*****************************************************************
 *
 *  Random simulation of the labels for Int_trees class
 *
 **/

void Int_trees::iid_simulation(const Distribution& distrib)
{

   int t, i; // var,
   vertex_iterator it, end;
   value v(_nb_integral, 0);

   if (trees != NULL)
   {
       for(t= 0; t < _nb_trees; t++)
       {
          Tree_tie::tie(it, end)= trees[t]->vertices();
          while ( it != end )
          {
             for(i= 0; i < _nb_integral; i++)
                    v.Int(i)= distrib.simulation();

             trees[t]->put(*it, v);
             it++;
          }

       }
       build_characteristics();
       build_size_frequency_distribution();
       build_nb_children_frequency_distribution();
   }
}

/*****************************************************************
 *
 *  Copy operator of Int_trees class
 *
 **/

// void Int_trees::copy(const Int_trees &otrees)
