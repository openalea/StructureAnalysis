/* -*-c++-*- 
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture 
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr) 
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


#ifndef SB_MATCHING_DISTANCE_TABLE
#define SB_MATCHING_DISTANCE_TABLE

#include "definitions.h"
#include "nodecost.h"
#include "mnodecost.h"
#include "wnodecost.h"
#include "tlnodecost.h"
#include "wnodecostnonconst.h"
#include "distancetable.h"
//#include <STAT/general.h>
//#include <STAT/vectors.h>

/**
 *\class MatchingDistanceTable
 *\brief Retains the distance between all the subtrees of two tre-graphs.
 * The distance is computing from a given NodeCost.
 *\par Requirements
 *  - Two TreeGraph;
 *  - A NodeCost
 *\author Pascal ferraro
 *\date 1999
 */

class MatchingDistanceTable
{
  friend class Matching;
  friend class MatchPath;

  public :
  MatchingDistanceTable(){};
  MatchingDistanceTable(TreeGraph& ,TreeGraph& ,NodeCost&);
  ~MatchingDistanceTable(){};
  void make(TreeGraph& ,TreeGraph& ,NodeCost& );
  DistanceType getDBT(int ,int ) const;
  DistanceType getDBF(int ,int ) const;
  DistanceType getDBOrderedF(int ,int ) ;
  void putDBT(int ,int ,DistanceType );
  void putDBF(int ,int ,DistanceType );
  void openDistancesVector(int );
  void closeDistancesVector(int );
  //Initialisation
  void putInputTreeToEmpty(int,DistanceType);
  DistanceType getInputTreeToEmpty(int);
  void putReferenceTreeFromEmpty(int,DistanceType);
  DistanceType getReferenceTreeFromEmpty(int);
  
  DistanceType inputTreeToEmpty(int );
  DistanceType inputForestToEmpty(int );
  DistanceType inputOrderedForestToEmpty(int );
  DistanceType referenceTreeFromEmpty(int );
  DistanceType referenceForestFromEmpty(int );
  DistanceType referenceOrderedForestFromEmpty(int );
  DistanceType getICost(int& ) const;
  DistanceType getDCost(int& ) const;
  DistanceType getCCost(int ,int ) const;
  NodeCostType getCostType() { return(ND->type()); };
  DistanceTable* getDistanceTable() {return(&_treeDistTable);};

  private :
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

