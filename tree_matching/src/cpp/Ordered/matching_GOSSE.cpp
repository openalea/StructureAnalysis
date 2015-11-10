#include "matching_GOSSE.h"


// -------------
// Constructeur
// -------------
MatchingGosse::MatchingGosse(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  MatchingGomse(input,reference,nodeDistance,0)
{
  TreeNode* node =  (&input)->getNode(0);
  int scale = (&input)->getMTG()->vscale(node->getVertex());
  _mss = new MS_O_Similarity(input, reference ,nodeDistance,scale);
}
