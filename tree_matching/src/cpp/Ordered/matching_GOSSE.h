//-*-c++-*-
#ifndef _MATCHING_GOSSE_H_
#define _MATCHING_GOSSE_H_
#include "matching_GOMSE.h"


class MatchingGosse : public MatchingGomse
{
public:
  MatchingGosse(TreeGraph&, TreeGraph&, NodeCost&);
};
#endif
