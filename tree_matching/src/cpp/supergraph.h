/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.Ferraro (pascal.ferraro@cirad.fr)
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



#ifndef SB_SUPERGRAPH_HEADER
#define SB_SUPERGRAPH_HEADER

#include "treegraph.h"

class SuperGraph {

public :
  SuperGraph(TreeGraph &, TreeGraph &);
  ~SuperGraph();
  int outdeg(node);
  std::list<node> getOutNodes(node);
  int getNbVertex(int ); // parameter indicates T1 or T2 (1 =>T1 and 2=>T2)
  int getVertexId(int);
 TreeGraph* T1;
  TreeGraph* T2;
private :
  // Id of nodes in T1 and T2 must be defined as follow :
  // Id of T1 nodes : from 0 to _size1 -1
  // Id of T2 nodes : from _size1 to _size1+_size2-1

 
  int _size1; // Size of T1
  int _size2; // Size of T2
};


#endif

