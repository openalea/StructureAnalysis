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


#ifndef SB_TREE_MATCH_HEADER
#define SB_TREE_MATCH_HEADER

#include "definitions.h"
#include "matching.h"
#include "self_similar_matching.h"
#include "end_space_free_matching.h"
#include "local_matching.h"
#include "matching3.h"
#include "matching4.h"
#include "matching_with_complex.h"
#include "self_similar_complex_matching.h"
#include "multiscale_matching.h"
#include "JiangWangZhangMetric2.h"
#include "TichitFerraro.h"
#include "TichitFerraro1.h"
#include "FerraroOuangraoua.h"
#include "FerraroOuangraoua1.h"
#include "FerraroOuangraoua2.h"
#include "matching_GPOSDE.h"
#include "selkow.h"
#include "matching2.h"
#include "matching_sequence.h"
//#include "extmatching.h"
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
//#include "STAT_TOOL/vectors.h"
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

enum MatchingType { TEST=0,
                    BY_WEIGHTS,
                    BY_TOPOLOGY,
                    BY_INCLUSION,
                    BY_COMPONENTS,
                    BY_EUCLIDIAN_DISTANCE,
                    BY_COMPLEX,
                    ORDERED_MATCHING,
		    SEQUENCE_MATCHING,
		    SELKOW_MATCHING,
		    END_SPACE_FREE,
		    FRACTAL,
		    EDITION,
                    ALIGNMENT,
                    SMALLESTCOMMONSUPERTREE,
                    LARGESTCOMMONSUBTREE

                  };

enum MappingType{ GLOBAL,
		  LOCAL
};

enum Mapping{ GENERAL,
	      MINIMUMCONNECTEDCOMPONENTS,
	      ENDSPACEFREE,
	      LCA,
              DISTANCE,
	      SIMILARITY
};

enum ScaleType{SINGLESCALE,
	       MULTISCALE,
	       QUOTIENTED
};

enum OrderedType { JIANG_WANG_ZHANG,
		   TICHIT_FERRARO,
		   TICHIT_FERRARO1,
		   FERRARO_OUANGRAOUA,
		   FERRARO_OUANGRAOUA1,
		   FERRARO_OUANGRAOUA2,
		   PARTIAL_ORDER,
		   ZHANG_SHASHA
                 };

enum InDelCostType {CONST_INDEL_COST,
                    NON_CONST_INDEL_COST};


/**
 *\class TreeMatch
 *\brief General class for computing a mapping between two tree graphs.
 *\author Pascal ferraro
 *\date 1999
 */


class TreeMatch
{
  public :

  TreeMatch();

  TreeMatch( MTG& mtg, Array* roots);

  TreeMatch( MTG& mtg,
             Array* roots,
             Array* local_functions,
             AMString matching_type,
             AMString ordered_type,
             int self_similarity,
             const Vector_distance &ivect,
             double coeff);

  TreeMatch( MTG& mtg,
             Array* roots,
             Array* local_functions,
             AMString matching_type,
             AMString ordered_type,
             int self_similarity,
             const Vector_distance &ivect,
             double coeff,NodeCost& Nd,int scale);

  ~TreeMatch();

  void weigthedMatching();

  void topologicalMatching();


  void localWeigthedMatching();

  void localTopologicalMatching();

  void selfSimilarityWeigthedMatching();

  void selfSimilarityTopologicalMatching();

  void selfSimilarityComplexMatching();

  void endSpaceFreeTopologicalMatching();
  void fractalTopologicalMatching();

   void sequenceMatching();

   void selkowMatching();

  //  void testMatching();


  DistanceType MatchByAttribute( TreeGraph& T1,
                                 TreeGraph& T2,
                                 Sequence * S);

  DistanceType LocalMatchByAttribute( TreeGraph& T1,
                                 TreeGraph& T2,
                                 Sequence * S);

  DistanceType MatchByAttribute( TreeGraph& T1,
                                 TreeGraph& T2);

  DistanceType JiangWangZhangMatching( TreeGraph& T1,
                                 TreeGraph& T2,
                                Sequence * S);

  DistanceType TichitFerraroMatching( TreeGraph& T1,
				      TreeGraph& T2,
				      Sequence * S,
				      int scale);

  DistanceType TichitFerraro1Matching( TreeGraph& T1,
				      TreeGraph& T2,
				      Sequence * S,
				      int scale);

  DistanceType FerraroOuangraouaMatching( TreeGraph& T1,
				      TreeGraph& T2,
				      Sequence * S,
				      int scale);

  DistanceType FerraroOuangraoua1Matching( TreeGraph& T1,
				      TreeGraph& T2,
				      Sequence * S,
				      int scale);

  DistanceType FerraroOuangraoua2Matching( TreeGraph& T1,
				      TreeGraph& T2,
				      Sequence * S,
				      int scale);

  DistanceType MatchPartialOrder( TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchTest( TreeGraph& T1,
                          TreeGraph& T2,
                          Sequence * S);


  DistanceType MatchByTopology( TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType LocalMatchByTopology(    TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchByTopology( TreeGraph& T1,
                                TreeGraph& T2);

  DistanceType MatchEndSpaceFreeByTopology( TreeGraph& T1,
					    TreeGraph& T2);

  DistanceType MatchByComplex ( TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchByComplex ( TreeGraph& T1,
                                TreeGraph& T2);

  DistanceType MatchByInclusion(TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchByComponents(TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchSequences( TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  DistanceType MatchSelkow( TreeGraph& T1,
                                TreeGraph& T2,
                                Sequence * S);

  std::ostream& introduce(std::ostream& out_info) const;

  std::ostream& viewAllMatching(std::ostream& out_fich) const;

  std::ostream& viewOneMatching(std::ostream& out_info,int imp_tree,int ref_tree) const;

  ostream& saveDistanceMatrix(ostream& out_fich) const;

  DistanceType viewDistanceMatching(std::ostream& out_info,int imp_tree,int ref_tree) const;
  DistanceType viewNormalizedDistance(std::ostream& out_info,int imp_tree,int ref_tree);
  bool mtg_write( const char *path ) const{
	  return _trees[0]->mtg_write(path);
  }

  void viewOrder();

  SLArray* getList(int i_tree,int r_tree);

  DistanceType  getDist(int i_tree,int r_tree);
  DistanceType  getNormalizedDistance(int i_tree,int r_tree);
  DistanceType  getSelfSimilarNormalizedDistance(int i_tree,int r_tree);
  int  viewVertex(int i_tree);


  void upDate(MTG&);

  TreeGraph* getTree(int ) const;

  Distance_matrix* getMatrix();

  Distance_matrix* getSelfSimilarityMatrix();

  DistanceType getDistance(int,int) const ;
  DistanceType getTime(int,int) const {
	  return 0;
  }

  DistanceType getSelfSimilarDistance(int,int) const{
	  return 0;
  }


  Sequence* getSequence(int,int) const ;

  Sequence*** getSequenceTab(int,int) const ;

  int getNbTree() const;

  MatchingType getMatchingType() const;

  NodeCostType getLocalDist() const;

  int getMaxOrder() const;

  int getNbVertexMax() const;
  DistanceType getParacladialCoefficient(int i_tree);
  SLArray* getParacladialCoefficients();
  int getParacladium(int i_tree);
  int getParacladium(int i_tree,int r_tree) const{
	  return 0;
  }
 SLArray* getParacladia();
 DistanceType getC(std::ostream& out_info,int imp_tree,int ref_tree) const{
 return 0;
 }
 
  DistanceType getBMCoefficient(int i_tree);
  SLArray* getBranchMappingCoefficients();
  int getBranchMapped(int i_tree);
  SLArray* getBranchMapped();

  DistanceType getBMNormalizedCoefficient(int i_tree);
  SLArray* getBranchMappingNormalizedCoefficients();
  int getNormalizedBranchMapped(int i_tree);
  SLArray* getNormalizedBranchMapped();






  protected:

  // data members;
    MTG*                _mtg;
    VIdList*            _roots;
    TreeVector          _trees;
    NodeFunctionList*   _localFun;
    MatchingType        _matchingType;
    OrderedType         _orderedType;
    int                 _selfSimilarity;
    DistanceVectorTable _distances;
    IntVectorTable      _intDistances;
    SequenceArray       _sequences;
    SequenceMatrix      _sequences_tab;
    int                 _maxOrder;
    int                 _maxNbVertex;
    char*               _fileName;
    int                 _nbTree;
    DistanceType        _InsDelCostCoeff;
    Vector_distance     _vectorDist;
    ValueVector         _dispersion;
    ValueVector         _maxValue;
    ValueVector         _minValue;
  NodeList            _trunk;
  DistanceVectorTable	_time;
  NodeCost* _nd;
  int _scale;
public:
  void putDistance(DistanceType,int,int);
  void putTime(DistanceType,int,int){
  }
 void putInt(int,int,int);
  void putDistanceMatrix(MatchingDistanceTable);
  void putSequence(Sequence*,int,int);
  void putSequenceTab(Sequence***,int,int);

  void standardization();
  DistanceType meanAbsoluteDifferenceComputation(int);


 };

#endif


