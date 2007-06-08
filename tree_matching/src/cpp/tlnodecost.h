
#ifndef SB_TL_NODE_COST_HEADER
#define SB_TL_NODE_COST_HEADER

#include "definitions.h"
#include "treenode.h"
#include "nodecost.h"


class TopologicalLocalNodeCost : public NodeCost
{
  public :

    /** Default constructor. */
    TopologicalLocalNodeCost(){}

    /** Destructor. */
    ~TopologicalLocalNodeCost(){}

    /** Constructs a NodeCost with the type /e type and with a Vector Distance, a dispersion, a maximum and minimum values.and a product coefficient for indel cost. */
    TopologicalLocalNodeCost(NodeCostType);


    /** Returns the insertion cost of /e node */
    DistanceType getInsertionCost(TreeNode* );

    /** Returns the deletion cost of /e node */
    DistanceType getDeletionCost(TreeNode* );

    /** Returns the changing cost between /e i_node and /e r_node*/
    DistanceType getChangingCost(TreeNode* ,TreeNode* );

    /** Returns the type of the node cost*/
    NodeCostType type() const { return(_type); }

    
};

#endif

