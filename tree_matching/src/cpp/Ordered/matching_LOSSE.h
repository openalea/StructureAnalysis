//-*-c++-*-
#ifndef _MATCHING_LOSSE_H_
#define _MATCHING_LOSSE_H_
#include "matching_LOMSE.h"


class MatchingLosse : public MatchingLomse
{
public:
  MatchingLosse(TreeGraph&, TreeGraph&, NodeCost&);
};
#endif
