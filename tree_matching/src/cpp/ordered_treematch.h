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


#ifndef SB_ORDERED_TREE_MATCH_HEADER
#define SB_ORDERED_TREE_MATCH_HEADER

#include "definitions.h"

/* Déclaration des differents matching */
#include "Ordered/matching_GOMDE.h"
#include "Ordered/matching_GOQDE.h"
#include "Ordered/matching_GOMSE.h"
#include "Ordered/matching_GOSDE_lca.h"
#include "Ordered/matching_LOMSE.h"
#include "Ordered/matching_LOSSE.h"
#include "Ordered/matching_GOSSE.h"
#include "Ordered/matching_O.h"


#include "matching.h"
#include "sequence.h"
#include "wnodecost.h"
#include "wnodecostnonconst.h"
#include "mnodecost.h"
#include "tlnodecost.h"
#include "treenode.h"
#include "amlnodefunction.h"
#include "mdtable.h"

#include "aml/array.h"
#include "aml/amstring.h"
#include "stat_tool/distance_matrix.h"

#include <iostream>
#include <time.h>

typedef std::vector<Sequence*>               SequenceVector;
typedef std::vector<SequenceVector*>         SequenceArray;
typedef std::vector<Sequence***>             SequenceTabVector;
typedef std::vector<SequenceTabVector*>      SequenceMatrix;
typedef std::vector<VId>                     VertexVector;
typedef std::list<VId>                       VertexList;
typedef std::list<int>                       IndexList;
typedef std::list<int>                       Component;
typedef std::list<Component*>                ComponentTable;
typedef std::vector<TreeGraph*>              TreeVector;

#include "treematch.h"
// enum MatchingType { EDITION,
//                     ALIGNMENT,
//                     SMALLESTCOMMONSUPERTREE,
//                     LARGESTCOMMONSUBTREE
// };






/**
 *\class TreeMatch
 *\brief General class for computing a mapping between two tree graphs.
 *\author Pascal ferraro
 *\date 1999
 */


class TreeMatch_O : public TreeMatch
{
  public :

  TreeMatch_O();

  TreeMatch_O( MTG& mtg, Array* roots);

  TreeMatch_O( MTG& mtg,
	       Array* roots,
	       Array* local_functions,
	       AMString matching_type,
	       AMString mapping_type,
	       AMString mapping,
	       AMString scale_type,
	       const stat_tool::VectorDistance &ivect,
	       double coeff);
    

  ~TreeMatch_O();

  DistanceType MatchByTopology(TreeGraph& Tree1,
			       TreeGraph& Tree2,
			       Sequence * S);
  void topologicalMatching();
  void weightedMatching();

  //Fonctions de ordered_matching_extract

//   stat_tool::DistanceMatrix* getMatrix();
//   SLArray* getList(int i_tree,int r_tree);
//   DistanceType getDist(int i_tree,int r_tree);
//   DistanceType getTime(int inp_tree,int ref_tree) const;
//   DistanceType getNormalizedDistance(int i_tree,int r_tree);
//   int viewVertex(int i_tree);
//   void putDistance(DistanceType distance,int inp_tree,int ref_tree);
//   void putTime(DistanceType time,int inp_tree,int ref_tree);
//   void putInt(int distance,int inp_tree,int ref_tree);
//   void putDistanceMatrix(MatchingDistanceTable matrix);
//   DistanceType getDistance(int inp_tree,int ref_tree) const;
//   void putSequence(Sequence* sequence,int inp_tree,int ref_tree);
//   void putSequenceTab(Sequence*** sequence_tab,int inp_tree,int ref_tree);
//   Sequence* getSequence(int inp_tree,int ref_tree) const;
//   Sequence*** getSequenceTab(int inp_tree,int ref_tree) const;
//   TreeGraph* getTree(int tree) const;
//   int  getNbTree() const;
//   MatchingType getMatchingType() const;
//   int  getMaxOrder() const;
//   int  getNbVertexMax() const;
//   bool mtg_write( const char *path ) const;

    private :


//   // data members;
//   MTG*                _mtg;
//   VIdList*            _roots;
//   TreeVector          _trees;
//   NodeFunctionList*   _localFun;
//   MatchingType        _matchingType;
   MappingType        _mappingType;
   Mapping            _mapping;
   ScaleType         _scaleType;
//   DistanceVectorTable _distances;
//   SequenceArray       _sequences;
//   SequenceMatrix      _sequences_tab;
//   int                 _maxOrder;
//   int                 _maxNbVertex;
//   char*               _fileName;
//   int                 _nbTree;
//   DistanceType        _InsDelCostCoeff;
//   stat_tool::VectorDistance _vectorDist;
//   ValueVector         _dispersion;
//   ValueVector         _maxValue;
//   ValueVector         _minValue;
//   NodeList            _trunk;
//   DistanceVectorTable	_time;
//   NodeCost* _nd;
//   int _scale;
//   IntVectorTable      _intDistances;
//   int                 _selfSimilarity;
//   // private functions
 
 


   void standardization();
   DistanceType meanAbsoluteDifferenceComputation(int);
 };

#endif


