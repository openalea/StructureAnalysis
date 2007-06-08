/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/matching_GPOSDE.h,v $
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


/*matching ne tenant pas compte du nombre de composantes connexes*/

#ifndef SB_MATCHING_GPOSDE_HEADER
#define SB_MATCHING_GPOSDE_HEADER
#include <vector>
#include "definitions.h"
#include "matchpath.h"
#include "choicetable.h"
#include "mdtable.h"
#include "treegraph.h"
#include "sequence.h"
#include "nodecost.h"
#include "wnodecost.h"
#include "mnodecost.h"


typedef std::vector<int> CaseVector;

/**
 *\class MatchingGposde With Complex
 *\brief Algorithm for comparing two unordered tree graph
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended
 * the algorithm proposed by Zhang \cite{Zha93}. The distance is computed as the minimum
 * cost of sequence of edit operations needed to transform one tree graph into another one.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal ferraro
 *\date 1999
 */

typedef  std::vector<std::vector<DistanceType > > DistanceTab;
typedef  std::vector<DistanceType> InDelTab;


class MatchingGposde
{
  
public :
  //Constructor
  MatchingGposde(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  
  ~MatchingGposde();
  
  //Distance Between Trees
  void distanceBetweenTree(int ,int ) ;
  //Distances Betweeen Forest Class
  void distanceBetweenForestClass(int ,int ) ;
  //Distances Betweeen Classes
  void distanceBetweenClass(int ,int ) ;
  
  

  DistanceType getInsertionCost(){return (getInBT(0,0));}
  DistanceType getDeletionCost(){return (getDeBT(0,0));}
  DistanceType getMatchingCost(){return (getSuBT(0,0));}
  Sequence* getSequence(){return _alignement;} 
  
  //Operator
  DistanceType getDBT(int ,int );
  DistanceType getDBFC(int ,int );
  DistanceType getDBC(int ,int );
  
  DistanceType getInsertCost(int);
  DistanceType getDeleteCost(int);
  DistanceType getSubstitutionCost(int,int);
  DistanceType  match();

  DistanceType getInBT(int ,int ) const;
  DistanceType getInBC(int ,int ) const;
  DistanceType getInBFC(int ,int ) const;
  DistanceType getDeBT(int ,int ) const;
  DistanceType getDeBC(int ,int ) const;
  DistanceType getDeBFC(int ,int ) const;
  DistanceType getSuBT(int ,int ) const;
  DistanceType getSuBC(int ,int ) const;
  DistanceType getSuBFC(int ,int ) const;
  int Lat(ChoiceList*, int );
  void TreeList(int ,int ,Sequence& );
  void ForestClassList(int ,int ,Sequence&);

protected :
  TreeGraph* T1;
  TreeGraph* T2;
  NodeCost* CostMatrix;
  Sequence* _alignement;
  
  DistanceTab _classes;
  DistanceTab _forestclasses;
  
  MatchingDistanceTable* _distances;
  InDelTab _insert;
  InDelTab _delete;
  
  MatchingDistanceTable* _insertCost;
  MatchingDistanceTable* _deleteCost;
  MatchingDistanceTable* _substitutionCost;
  MatchingDistanceTable* _insertClassCost;
  MatchingDistanceTable* _deleteClassCost;
  MatchingDistanceTable* _substitutionClassCost;
  ChoiceTable _choices;
  ChoiceTable _choicesFClass;
  int M(int,int);
  
};


#endif

