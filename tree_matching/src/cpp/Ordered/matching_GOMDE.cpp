#include "matching_GOMDE.h"
#include <iostream>
#include <cstdio>
#include <cmath>

using namespace std;

class Min {
  DistanceType value;
public:
  Min(DistanceType v):value(v) {}
  Min &operator<<(DistanceType v) {
    if (v<value) value=v;
    return *this;
  }
  operator DistanceType() {
    return value;
  }
};



// -------------
// Constructeur
// -------------

MatchingGomde::MatchingGomde(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance, int scale):
  Matching_O(input,reference,nodeDistance)
{

  _msm = new MS_O_Matching(input, reference ,nodeDistance,scale);

}

// -------------
// Destructeur
// -------------
MatchingGomde::~MatchingGomde()
{
  delete (MS_O_Matching*) _msm;

}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
/*DistanceType  MatchingGomde::run()
{

  return _msm->match();

}*/

/*DistanceType MatchingGomde::getInsertionCost(){return _msm->getInsertionCost();}
  DistanceType MatchingGomde::getDeletionCost(){return _msm->getDeletionCost();}
  DistanceType MatchingGomde::getMatchingCost(){return _msm->getMatchingCost();}

 Sequence* MatchingGomde::getSequence(){return _msm->getSequence();} 
*/
