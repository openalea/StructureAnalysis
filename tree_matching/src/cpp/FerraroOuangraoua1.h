//-*-c++-*-
#ifndef _FERRAROOUANGRAOUA1_H
#define _FERRAROOUANGRAOUA1_H

#include "MSSimilarity.h"


class FerraroOuangraoua1 {
public:
  FerraroOuangraoua1(TreeGraph&, TreeGraph&, NodeCost&, int scale);
  ~FerraroOuangraoua1();
  DistanceType run();
  int getvMax(){return 0;}
  int getwMax(){return 0;}
  float getInsertionCost();
  float getDeletionCost();
  float getMatchingCost();
  Sequence* getSequence(); 
private:
  MSSimilarity* _msm;
};

#endif



