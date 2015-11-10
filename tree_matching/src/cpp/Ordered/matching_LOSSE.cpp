#include "matching_LOSSE.h"


// -------------
// Constructeur
// -------------
MatchingLosse::MatchingLosse(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  MatchingLomse(input,reference,nodeDistance,0)
{
  TreeNode* node =  (&input)->getNode(0);
  int scale = (&input)->getMTG()->vscale(node->getVertex());
  _msl = new MSL_O_Similarity(input, reference ,nodeDistance,scale);
  
}
