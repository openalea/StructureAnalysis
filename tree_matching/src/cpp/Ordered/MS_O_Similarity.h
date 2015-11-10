//-*-c++-*-
#ifndef _MS_O_SIMILARITY_H
#define _MS_O_SIMILARITY_H

#include "MS_O_Matching.h"
#include "../treegraph.h" 
#include "../nodecost.h"
#include "../distancetable.h"
#include "../sequence.h"
#include <string>
using std::string;
#include <iostream>
#include <map>
#include <string>
using namespace std;



/**
 *\class MS Similarity
 *\To compare recursively Forests and subtrees of two quotiented tree graph.
 * In order to compute global similarities between Forest, we have extended 
 * the algorithm proposed by K.Zhang and A.Shapiro \cite{Zha93}. 
 * The distance is computed recursively using similarities yet computed.
 * We use two scales of decomposition.
 *\author Aïda Ouangraoua
 *\date 12/2004
 */


class MS_O_Similarity {
  
public:
 
  MS_O_Similarity (TreeGraph&, TreeGraph&, NodeCost&, int);
  ~MS_O_Similarity();
  DistanceType match();
  DistanceType getInsertionCost(){return insertCost;}
  DistanceType getDeletionCost(){return deleteCost;}
  DistanceType getMatchingCost(){return matchCost;}
  Sequence* getSequence(){return alignement;} 
  void printForest(int, int);  
  
  
  TreeGraph* T1;
  TreeGraph* T2;
  std::map<string, int> symbols;
  NodeCost* CostMatrix;
  static DistanceType CostMatrix1[NB_SYMBOLS][NB_SYMBOLS];
  
  int _init;
  
  std::vector<int> L1;
  std::vector<int> keyroots1;
  int sizeKeyroots1;
  std::vector<int> L2;
  std::vector<int> keyroots2;
  int sizeKeyroots2;
  
  void computeLeft(int numT, int root);
  void computeKeyroots(int numT);
  void precompute();
  
  ForestDistanceTable fDist;
  void computeForestsDistances(int , int); 
  int fix(int, int);
  int idx(int, int);

  DistanceType insLocalNode(TreeGraph*, int) ;
  DistanceType delLocalNode(TreeGraph*, int) ;
  DistanceType matchLocalNode(int, int) ;
  
  //  TreeDistanceTable _dst1;
  void computeTreeDist(int, int);

  DistanceType getTreeDistance(int, int);  
  
  Sequence* alignement;
  void makeAlignement(int, int);
  void subAlign(int, int, int, int);  
 
  
  DistanceType insertCost;
  DistanceType deleteCost;
  DistanceType matchCost; 
  int getKey1(int);
  int getKey2(int);
 
  static const DistanceType epsilon;
 
};

#endif
