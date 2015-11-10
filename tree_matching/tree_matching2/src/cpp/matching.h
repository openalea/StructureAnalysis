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
 *       $Id: matching.h 3258 2007-06-06 13:18:26Z dufourko $
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

#ifndef SB_MATCHING_HEADER
#define SB_MATCHING_HEADER

#include <vector>
#include "definitions.h"
#include "matchpath.h"
#include "choicetable.h"
#include "mdtable.h"
#include "treegraph.h"
#include "sequence.h"
#include "nodecost.h"



/**
 *\class Matching 
 *\brief Algorithm for comparing two unordered tree graph
 *\par Presentation
 * The distance between two trees is computed as the minimum
 * cost of sequence of edit operations needed to transform one tree graph into another one.
 *\par Requirements
 * - Two TreeGraphs defined at the same scale;
 * - A NodeCost (method for computing the local distance),
 *\author Pascal ferraro
 *\date 2009
 */



class TREEMATCH_API Matching
{

public :
  //Constructor
  Matching() {};

  Matching(const TreeGraphPtr& , const TreeGraphPtr& ,const NodeCostPtr&, MDTableType mdtable_type = STD ) ;

  //Destructor
  virtual ~Matching();
  
  //Distance Between Trees
  virtual DistanceType distanceBetweenTree(int ,int ) ;
  //Distances Betweeen Forest
  virtual DistanceType distanceBetweenForest(int ,int ) ;
  
  
  
  //Operator
  inline DistanceVectorTable getDistanceTable(){ 
    if (_distances->getType() == STD)
      return ((StdMatchingDistanceTable*)_distances)->getDistanceTable(); 
    else
      exit(0);
  }

  DistanceType getDBT(int ,int ) const;
  DistanceType getDBF(int ,int ) const;
  virtual DistanceType  match();

  void getList(int ,int ,Sequence*);
  // void ForestList(int ,int ,Sequence& );
  // void TreeList(int ,int ,Sequence& );

  int operator()(int );
  // int Lat(ChoiceList* L, int vertex);
  

  inline const ChoiceTable& getChoiceTable() const { return _choices; }

  bool verbose;

protected :
  TreeGraphPtr T1;
  TreeGraphPtr T2;
  MatchingDistanceTable * _distances;
  ChoiceTable _choices;
  NodeCostPtr ND;
  MatchPath _restrMapp;
  VertexVector _restrMappList;

  int M(int,int);

  MDTableType _mdtable_type;
  
};

class TREEMATCH_API ExtMatching : public Matching
{

public :
  ExtMatching(const TreeGraphPtr& input,const TreeGraphPtr& reference,const NodeCostPtr& nodeDistance, MDTableType mdtable_type  ):Matching(input,reference,nodeDistance,mdtable_type){};
  virtual ~ExtMatching();
  
  //Distance Between Trees
  virtual DistanceType distanceBetweenTree(int ,int ) ;
  
};

#endif

