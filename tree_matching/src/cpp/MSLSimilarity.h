//-*-c++-*-
#ifndef _MSLSIMILARITY_H
#define _MSLSIMILARITY_H

#include "MSSimilarity.h"
#include "treegraph.h" 
#include "nodecost.h"
#include "distancetable.h"
#include "sequence.h"
#include <string>
using std::string;
#include <iostream>
#include <map>
#include <string>
using namespace std;



/**
 *\class MSL Similarity
 *\To compare recursively Forests and subtrees of two quotiented tree graph.
 * In order to compute local similarities between Forest, we have extended 
 * the algorithm proposed by K.Zhang and A.Shapiro \cite{Zha93}. 
 * The distance is computed recursively using similarities yet computed.
 * We use two scales of decomposition.
 *\author Aïda Ouangraoua
 *\date 12/2004
 */




class MSLSimilarity {
  
public:
 
  MSLSimilarity (TreeGraph&, TreeGraph&, NodeCost&, int);
  ~MSLSimilarity();
  DistanceType match();
  float getInsertionCost(){return insertCost;}
  float getDeletionCost(){return deleteCost;}
  float getMatchingCost(){return matchCost;}
  Sequence* getSequence(){return alignement;} 
  int getvMax(){return vmax;}
  int getwMax(){return wmax;}
  
private: 
  DistanceType lsmax;
  int vmax;
  int wmax;
  int lvmax;
  int lwmax;
  TreeGraph* T1;
  TreeGraph* T2;
  std::map<string, int> symbols;
  NodeCost* CostMatrix;
  static float CostMatrix1[NB_SYMBOLS][NB_SYMBOLS];
  
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
  
  float insLocalNode(TreeGraph*, int) ;
  float delLocalNode(TreeGraph*, int) ;
  float matchLocalNode(int, int) ;
  
  TreeDistanceTable _dst1;
  void computeTreeDist(int, int);
  
  float getTreeDistance(int, int);  
  
  Sequence* alignement;
  void makeAlignement(int, int);
  void subAlign(int, int, int, int);  
  
  
  DistanceType insertCost;
  DistanceType deleteCost;
  DistanceType matchCost;
  int getKey1(int);
  int getKey2(int);
  
  static const float epsilon;
};

#endif
