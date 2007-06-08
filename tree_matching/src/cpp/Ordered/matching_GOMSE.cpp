#include "matching_GOMSE.h"
#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;




// -------------
// Constructeur
// -------------
MatchingGomse::MatchingGomse(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale):
  Matching_O(input,reference,nodeDistance)
{ 
  _mss = new MS_O_Similarity(input, reference ,nodeDistance,scale); 
}

// -------------
// Destructeur
// -------------
MatchingGomse::~MatchingGomse()
{
  delete (MS_O_Similarity*) _mss;
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  MatchingGomse::run()
{

  return _mss->match();

}

DistanceType MatchingGomse::getInsertionCost(){return _mss->getInsertionCost();}
DistanceType MatchingGomse::getDeletionCost(){return _mss->getDeletionCost();}
DistanceType MatchingGomse::getMatchingCost(){return _mss->getMatchingCost();}

Sequence* MatchingGomse::getSequence(){return _mss->getSequence();} 
