/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
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



#ifndef SB_LOCALMATCHING_U_HEADER
#define SB_LOCALMATCHING_U_HEADER

#include <vector>
#include "../definitions.h"
#include "../localmatchpath.h"
#include "../choicetable.h"
#include "../mdtable.h"
#include "../treegraph.h"
#include "../sequence.h"
#include "../nodecost.h"
#include "../wnodecost.h"
#include "../mnodecost.h"
#include "matching_unordered.h"

typedef std::vector<int> CaseVector;


/**
 *\class Local Matching
 *\brief Algorithm for evaluating local similarity between two unordered tree graph
 *\par Presentation
 * In order to compute a distance between two multiscale tree graph, we have extended
 * the algorithm proposed by Zhang \cite{Zha93} and Smith and Waterman \cite{smi81}. The distance is computed as the minimum
 * cost of sequence of edit operations needed to transform one tree graph into another one.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal ferraro
 *\date 2003
 */

class LocalMatching_U :public Matching_U
{
public :
  //Constructor
  LocalMatching_U() {};
  LocalMatching_U(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  void make(TreeGraph& , TreeGraph& ,NodeCost& ) ;
  //Destructor
  ~LocalMatching_U();

  //Distance Between Trees
  virtual DistanceType distanceBetweenTree(int ,int ) ;
  //Distances Betweeen Forest
  virtual DistanceType distanceBetweenForest(int ,int ) ;



  //Operator
  DistanceType getDBT(int ,int ) const;
  DistanceType getDBF(int ,int ) const;
  DistanceType  match();
  DistanceType getDistanceMatrix(int ,int) const;
  void TreeList(int ,int ,Sequence& );
  int operator()(int );
  int getNbMin(){return _nb_min;};
  int getMaxRef();
  int getMaxInput();

protected :
  LocalMatchPath _localrestrMapp;
  double _maxScore;
  int _maxInput;
  int _maxRef;
  int _init;
};


#endif

