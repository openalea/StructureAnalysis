/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source$
 *       $Id: mdtable.h 3258 2007-06-06 13:18:26Z dufourko $
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


#ifndef SB_MATCHING_DISTANCE_TABLE
#define SB_MATCHING_DISTANCE_TABLE

#include "definitions.h"
#include "nodecost.h"
#include "treegraph.h"
/**
 *\class MatchingDistanceTable
 *\brief Retains the distance between all the subtrees of two tre-graphs.
 * The distance is computing from a given NodeCost.
 *\par Requirements
 *  - Two TreeGraph;
 *  - A NodeCost
 *\author Pascal ferraro
 *\date 2009
 */

typedef std::vector<double> DistanceVector;
typedef std::vector<DistanceVector> DistanceTable;

class MatchingDistanceTable
{
//   friend class Matching;
//   friend class MatchPath;

  public :
  MatchingDistanceTable(){};
  MatchingDistanceTable(TreeGraph& ,TreeGraph& ,NodeCost&);
  ~MatchingDistanceTable(){};
  DistanceType getDBT(int ,int ) const;
  DistanceType getDBF(int ,int ) const;
  void putDBT(int ,int ,DistanceType );
  void putDBF(int ,int ,DistanceType );
  //Initialisation
  void putInputTreeToEmpty(int,DistanceType);
  DistanceType getInputTreeToEmpty(int);
  void putReferenceTreeFromEmpty(int,DistanceType);
  DistanceType getReferenceTreeFromEmpty(int);
  
  DistanceType inputTreeToEmpty(int );
  DistanceType inputForestToEmpty(int );
  DistanceType referenceTreeFromEmpty(int );
  DistanceType referenceForestFromEmpty(int );
  DistanceType getICost(int& ) const;
  DistanceType getDCost(int& ) const;
  DistanceType getCCost(int ,int ) const;
  NodeCostType getCostType() { return(ND->type()); };
  DistanceTable &getDistanceTable() {return(_treeDistTable);};

  private :
  int _n1;
  int _n2;
  TreeGraph* T1;
  TreeGraph* T2;
  NodeCost* ND;

  protected :
  DistanceVector _inputTreeToEmpty;
  DistanceVector _referenceTreeFromEmpty;
  DistanceTable _treeDistTable;
  DistanceTable _forestDistTable;
};

#endif 

