//-*-c++-*-
#ifndef _FERRAROOUANGRAOUA_H
#define _FERRAROOUANGRAOUA_H

#include "MSMatching.h"


class FerraroOuangraoua {
public:
  FerraroOuangraoua(TreeGraph&, TreeGraph&, NodeCost&, int scale);
  ~FerraroOuangraoua();
  DistanceType run();
 int getvMax(){return 0;}
  int getwMax(){return 0;}
  float getInsertionCost();
  float getDeletionCost();
  float getMatchingCost();
  Sequence* getSequence(); 
  MSMatching* getMSM(){return _msm;};
private:
  MSMatching* _msm;
};

#endif



