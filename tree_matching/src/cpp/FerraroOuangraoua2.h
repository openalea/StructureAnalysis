//-*-c++-*-
#ifndef _FERRAROOUANGRAOUA2_H
#define _FERRAROOUANGRAOUA2_H

#include "MSLSimilarity.h"


class FerraroOuangraoua2 {
public:
  FerraroOuangraoua2(TreeGraph&, TreeGraph&, NodeCost&, int scale);
  ~FerraroOuangraoua2();
  DistanceType
  run();
  int getvMax(){return _msm->getvMax();}
  int getwMax(){return _msm->getwMax();}
  float getInsertionCost();
  float getDeletionCost();
  float getMatchingCost();
  Sequence* getSequence(); 
  
private:
  MSLSimilarity* _msm;
};

#endif



