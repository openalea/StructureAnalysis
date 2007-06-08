/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (aida.ouangraoua@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/Ordered/matching_GOQDE.h,v $
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

#ifndef SB_MATCHING_GOQDE_HEADER
#define SB_MATCHING_GOQDE_HEADER
#include <vector>
#include "../definitions.h"
#include "../choicetable.h"
#include "../mdtable.h"
#include "../treegraph.h"
#include "../sequence.h"
#include "../nodecost.h"
#include "../wnodecost.h"
#include "../mnodecost.h"
#include "matching_O.h"


typedef std::vector<int> CaseVector;

/**
 *\class MatchingGoqde
 *\brief Algorithm for comparing two quotiented ordered tree graph
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended
 * the algorithm proposed by Zhang \cite{Zha93}. The distance is computed as the minimum
 * cost of sequence of edit operations needed to transform one tree graph into another one.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Aida Ouangraoua
 *\date 2006
 */

typedef  std::vector<std::vector<DistanceType> > DistanceTab;
typedef  std::vector<DistanceType> InDelTab;


class MatchingGoqde : public Matching_O
{
  
public :
  //Constructor
  MatchingGoqde(TreeGraph& , TreeGraph& ,NodeCost& ) ;



  ~MatchingGoqde();
  
 
  //Distances Betweeen Forest
  void distanceBetweenForest(int ,int) ;
  
  void computeDist(int , int , int , int , int , int );
  
  //Operator
  
  DistanceType getDBT(int ,int, int , int );
  DistanceType getDist(int ,int, int , int, int, int);
  DistanceType getFromFDist(int ,int, int , int, int, int);
  
  DistanceType getInsertCost(int j){return CostMatrix->getInsertionCost(T2->getNode(j));}
  DistanceType getDeleteCost(int i ){return CostMatrix->getDeletionCost(T1->getNode(i));}
  DistanceType getSubstitutionCost(int i,int j){return CostMatrix->getChangingCost(T1->getNode(i),T2->getNode(j));}
  
  void computeLeft(int , int);
  
  void computeKeyroots(int);
  
  void precompute();
  
  
  int fix(int, int);
  
  int idx(int, int);
  
  //void makeAlignement(int, int);
  //void subAlign(int, int, int, int);  
  
  int getKey1(int);
  int getKey2(int);
  

  int lca(int, int, int);
  
  int classRoot(int, int);
  
  
  DistanceType match();
  
  
 protected:
  TreeGraph* T1;
  TreeGraph* T2;
  NodeCost* CostMatrix;
  
  
  std::vector<int> L1;
  
  std::vector<int> L2;
  
  
  std::vector<int> keyroots1;
  std::vector<int> keyclasses1;
  int sizeKeyroots1;
  std::vector<int> keyroots2;
  std::vector<int> keyclasses2;
  int sizeKeyroots2;
  
  
  ForestDistanceTable fDist;
  ForestDistanceTable fDistC1;
  ForestDistanceTable fDistC2;
  ForestDistanceTable fDistC1C2;
  
  
  
  Sequence* alignement;
  
  
  int insertCost;
  int deleteCost;
  int matchCost; 
  
  static const float epsilon;
};


#endif

