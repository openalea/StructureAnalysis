#include "FerraroOuangraoua.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>

using namespace std;

class Min {
  float value;
public:
  Min(float v):value(v) {}
  Min &operator<<(float v) {
    if (v<value) value=v;
    return *this;
  }
  operator float() {
    return value;
  }
};



// -------------
// Constructeur
// -------------
FerraroOuangraoua::FerraroOuangraoua(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale)
{

  _msm = new MSMatching(input, reference ,nodeDistance,scale);

}

// -------------
// Destructeur
// -------------
FerraroOuangraoua::~FerraroOuangraoua()
{
  delete (MSMatching*) _msm;
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  FerraroOuangraoua::run()
{

  return _msm->match();

}

 float FerraroOuangraoua::getInsertionCost(){return _msm->getInsertionCost();}
  float FerraroOuangraoua::getDeletionCost(){return _msm->getDeletionCost();}
  float FerraroOuangraoua::getMatchingCost(){return _msm->getMatchingCost();}

 Sequence* FerraroOuangraoua::getSequence(){return _msm->getSequence();} 
