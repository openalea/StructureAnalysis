#include "FerraroOuangraoua2.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>

using namespace std;




// -------------
// Constructeur
// -------------
FerraroOuangraoua2::FerraroOuangraoua2(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale)
{
 _msm = new MSLSimilarity(input, reference ,nodeDistance,scale);
  
}

// -------------
// Destructeur
// -------------
FerraroOuangraoua2::~FerraroOuangraoua2()
{
  delete (MSLSimilarity*) _msm;
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------

DistanceType  FerraroOuangraoua2::run()
{
  
  return _msm->match();
  
}

float FerraroOuangraoua2::getInsertionCost(){return _msm->getInsertionCost();}
float FerraroOuangraoua2::getDeletionCost(){return _msm->getDeletionCost();}
float FerraroOuangraoua2::getMatchingCost(){return _msm->getMatchingCost();}

Sequence* FerraroOuangraoua2::getSequence(){return _msm->getSequence();} 

