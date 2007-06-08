//-*-c++-*-
#ifndef _MATCHING_GOMSE_H
#define _MATCHING_GOMSE_H

#include "matching_O.h"
#include "MS_O_Similarity.h"


class MatchingGomse : public Matching_O
{
public:
  MatchingGomse(TreeGraph&, TreeGraph&, NodeCost&, int);
  ~MatchingGomse();
  DistanceType run();
  DistanceType match(){return run();}
  int getvMax(){return 0;}
  int getwMax(){return 0;}
  DistanceType getInsertionCost();
  DistanceType getDeletionCost();
  DistanceType getMatchingCost();
  Sequence* getSequence(); 
  //void TreeList(int i,int j,Sequence& s){s = *(getSequence());}
protected:
  MS_O_Similarity* _mss;
};

#endif



