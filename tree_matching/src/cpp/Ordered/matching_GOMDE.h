//-*-c++-*-
#ifndef _MATCHING_GOMDE_H
#define _MATCHING_GOMDE_H
#include "matching_O.h"


class MatchingGomde : public Matching_O
{
public:
  MatchingGomde(TreeGraph&, TreeGraph&, NodeCost&, int scale);
  ~MatchingGomde();
  //DistanceType run();
  //DistanceType match(){return run();}
  int getvMax(){return 0;}
  int getwMax(){return 0;}
  //float getInsertionCost();
  //float getDeletionCost();
  //float getMatchingCost();
  //Sequence* getSequence(); 
  //void TreeList(int i,int j,Sequence& s){s = *(getSequence());}  
};

#endif



