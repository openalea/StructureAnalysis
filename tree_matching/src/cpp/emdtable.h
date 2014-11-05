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


#ifndef SB_EXTENDED_MATCHING_DISTANCE_TABLE
#define SB_EXTENDED_MATCHING_DISTANCE_TABLE

#include "definitions.h"
#include "nodecost.h"
#include "mnodecost.h"
#include "wnodecost.h"
#include "distancetable.h"


class ExtendedMatchingDistanceTable
{
  friend class ExtendedMatching;
  friend class MatchPath;
  
  public :

  ExtendedMatchingDistanceTable(){};
  ExtendedMatchingDistanceTable(TreeGraph& , TreeGraph& ,NodeCost&);
  ~ExtendedMatchingDistanceTable(){};
  void make(TreeGraph& , TreeGraph& ,NodeCost& );
  DistanceType getDBT(int ,int ) const ;
  DistanceType getDBF(int ,int ) const ;
  DistanceType getDBMF(int ,int ) const;
  DistanceType getDBPF(int ,int ) const;
  DistanceType getDBIMFRF(int ,int ) const;
  DistanceType getDBIFRMF(int ,int ) const;
  void putDBT(int ,int ,DistanceType );
  void putDBF(int ,int ,DistanceType );
  void putDBMF(int ,int ,DistanceType );
  void putDBPF(int ,int ,DistanceType );
  void putDBIMFRF(int ,int ,DistanceType  );
  void putDBIFRMF(int ,int  ,DistanceType );
  void openDistancesVector(int );
  void closeDistancesVector(int );

  //Initialisation
  DistanceType inputTreeToEmpty(int );
  DistanceType inputForestToEmpty(int );
  DistanceType referenceTreeFromEmpty(int );
  DistanceType referenceForestFromEmpty(int );
  DistanceType getICost(int ) const;
  DistanceType getDCost(int ) const;
  DistanceType getCCost(int ,int ) const;
  NodeCostType getCostType() { return(ND->type()); };
  DistanceTable* getDistanceTable() { return(&_treeDistTable); };
  private :

  TreeGraph* T1;
  TreeGraph* T2;
  NodeCost* ND;

  protected :
  DistanceVector _inputTreeToEmpty;
  DistanceVector _referenceTreeFromEmpty;
  DistanceTable _treeDistTable;
  DistanceTable _forestDistTable;
  DistanceTable _minusForestDistTable;
  DistanceTable _plusForestDistTable;
  DistanceTable _imfrfDistTable;
  DistanceTable _ifrmfDistTable;
  AmlBoolean sameComplex(TreeNode* ,TreeNode* ) const;

};

#endif 




