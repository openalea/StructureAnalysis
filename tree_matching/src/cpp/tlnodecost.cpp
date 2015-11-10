

#include"tlnodecost.h"

TopologicalLocalNodeCost::TopologicalLocalNodeCost( NodeCostType type)
{
  _type=type;
  _norm = L1;
}


DistanceType TopologicalLocalNodeCost::getInsertionCost(TreeNode* node)
{
  DistanceType cost;

  if (_norm == L1)
   {
     cost=node->getValue();
   }
  return(-cost);
}

DistanceType TopologicalLocalNodeCost::getDeletionCost(TreeNode* node)
{
  DistanceType cost;
  if(_norm == L1)
   {
     cost=node->getValue();
   }
  return(-cost);
}

DistanceType TopologicalLocalNodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost;
  if(_norm == L1)
   {
     //cost=i_node->getValue();
 cost=2-ABS(i_node->getValue()-r_node->getValue());
  }
  return(cost);
}


