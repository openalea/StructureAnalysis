#include "matching_LOMSE.h"
#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;




// -------------
// Constructeur
// -------------
MatchingLomse::MatchingLomse(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale):
  Matching_O(input,reference,nodeDistance)
{
  _msl = new MSL_O_Similarity(input, reference ,nodeDistance,scale);  
}

// -------------
// Destructeur
// -------------
MatchingLomse::~MatchingLomse()
{
  delete (MSL_O_Similarity*) _msl;

}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------

DistanceType  MatchingLomse::run()
{
  
  return _msl->match();
  
}

DistanceType MatchingLomse::getInsertionCost(){return _msl->getInsertionCost();}
DistanceType MatchingLomse::getDeletionCost(){return _msl->getDeletionCost();}
DistanceType MatchingLomse::getMatchingCost(){return _msl->getMatchingCost();}

Sequence* MatchingLomse::getSequence(){return _msl->getSequence();} 

