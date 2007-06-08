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


#ifndef SB_MATCHING_INT_TABLE
#define SB_MATCHING_INT_TABLE

#include "definitions.h"
#include "nodecost.h"
#include "inttable.h"
#include "treegraph.h"
//#include <STAT/general.h>
//#include <STAT/vectors.h>

/**
 *\class MatchingIntTable
 *\brief Retains the distance between all the subtrees of two tre-graphs.
 * The distance is computing from a given NodeCost.
 *\par Requirements
 *  - Two TreeGraph;
 *  - A NodeCost
 *\author Pascal ferraro
 *\date 2003
 */

class MatchingIntTable
{
  friend class Matching;
  friend class MatchPath;

  public :
  MatchingIntTable(){};
  MatchingIntTable(TreeGraph& ,TreeGraph& ,NodeCost&);
  ~MatchingIntTable(){};
  void make(TreeGraph& ,TreeGraph& ,NodeCost& );
  IntType getDBT(int ,int ) const ;
  IntType getDBF(int ,int ) const ;
  void putDBT(int ,int ,IntType );
  void putDBF(int ,int ,IntType );
  void openIntVector(int );
  void closeIntVector(int );
  //Initialisation
  IntType inputTreeToEmpty(int );
  IntType inputForestToEmpty(int );
  IntType referenceTreeFromEmpty(int );
  IntType referenceForestFromEmpty(int );
  IntType getICost(int& ) const;
  IntType getDCost(int& ) const;
  IntType getCCost(int ,int ) const;
  NodeCostType getCostType() { return(ND->type()); }
  IntTable* getIntTable() {return(&_treeDistTable);}

  private :
  TreeGraph* T1;
  TreeGraph* T2;
  NodeCost* ND;

  protected :
  IntVector _inputTreeToEmpty;
  IntVector _referenceTreeFromEmpty;
  IntTable _treeDistTable;
  IntTable _forestDistTable;
};

#endif 

