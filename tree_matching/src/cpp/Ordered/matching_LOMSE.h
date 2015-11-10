//-*-c++-*-
#ifndef _MATCHING_LOMSE_H_
#define _MATCHING_LOMSE_H_
#include "matching_O.h"
#include "MSL_O_Similarity.h"


class MatchingLomse : public Matching_O
{
public:
  MatchingLomse(TreeGraph&, TreeGraph&, NodeCost&, int scale);
  ~MatchingLomse();
  DistanceType run();
  DistanceType match(){return run();}
  int getvMax(){return _msl->getvMax();}
  int getwMax(){return _msl->getwMax();}
  DistanceType getInsertionCost();
  DistanceType getDeletionCost();
  DistanceType getMatchingCost();
  Sequence* getSequence(); 
  // void TreeList(int i,int j,Sequence& s){s = *(getSequence());}  
protected:
  MSL_O_Similarity* _msl;
};

#endif



