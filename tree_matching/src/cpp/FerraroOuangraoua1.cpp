#include "FerraroOuangraoua1.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>

using namespace std;




// -------------
// Constructeur
// -------------
FerraroOuangraoua1::FerraroOuangraoua1(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale)
{

  _msm = new MSSimilarity(input, reference ,nodeDistance,scale);

}

// -------------
// Destructeur
// -------------
FerraroOuangraoua1::~FerraroOuangraoua1()
{
  delete (MSSimilarity*) _msm;
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  FerraroOuangraoua1::run()
{

  return _msm->match();

}

 float FerraroOuangraoua1::getInsertionCost(){return _msm->getInsertionCost();}
  float FerraroOuangraoua1::getDeletionCost(){return _msm->getDeletionCost();}
  float FerraroOuangraoua1::getMatchingCost(){return _msm->getMatchingCost();}

 Sequence* FerraroOuangraoua1::getSequence(){return _msm->getSequence();} 
