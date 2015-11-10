//-*-c++-*-
#ifndef _MS_O_MATCHING_H
#define _MS_O_MATCHING_H

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

#define NB_SYMBOLS 24

/**
 *\class MS Matching
 *\To compare recursively Forests and subtrees of two quotiented tree graph.
 * In order to compute distances between Forest, we have extended 
 * the algorithm proposed by K.Zhang and A.Shapiro \cite{Zha93}. 
 * The distance is computed recursively using distances yet computed.
 * We use two scales of decomposition.
 *\author Aïda Ouangraoua
 *\date 04/2004
 */

typedef  std::vector<std::vector<std::vector<DistanceType> > > ForestDistanceTable;
typedef  std::vector<std::vector<DistanceType> > TreeDistanceTable;

class MS_O_Matching {
  
public:
 
  MS_O_Matching (TreeGraph&, TreeGraph&, NodeCost&, int);
  ~MS_O_Matching();
  DistanceType match();
  DistanceType getInsertionCost(){return insertCost;}
  DistanceType getDeletionCost(){return deleteCost;}
  DistanceType getMatchingCost(){return matchCost;}
  Sequence* getSequence(){return alignement;} 
  
  
  

  TreeGraph* T1;
  TreeGraph* T2;
  std::map<string, int> symbols;
  std::map<string, int> symbols2;
  NodeCost* CostMatrix;
  static DistanceType CostMatrix1[NB_SYMBOLS][NB_SYMBOLS];
  static DistanceType CostMatrix2[8][8];
  
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

  int getNbVertexAtScaleInit(TreeGraph*, int) ;
  
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
