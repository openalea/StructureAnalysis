from openalea.tree_matching import *

class MyNodeCost (NodeCost):
  def __init__(self):
     NodeCost.__init__(self)
  def getInsertionCost(self,a):
    print('ic')
    print(a)
    return 18
  def getDeletionCost(self,a):
    print('dc')
    print(a)
    return 50


def test_matching():
    tree1 = TreeGraph()
    tree1.addNode(0,-1)
    tree1.addNode(1,0)
    tree1.addNode(2,1)
    tree1.addNode(3,1)
    tree1.addNode(4,0)
    
    tree2 = TreeGraph()
    tree2.addNode(0,-1)
    tree2.addNode(1,0)
    tree2.addNode(2,0)
    tree2.addNode(3,0)
    
    node_cost = MyNodeCost()
    print(node_cost.getInsertionCost(tree1.getNode(4)))
    print(node_cost.getDeletionCost(tree1.getNode(4)))
    print('matching')
    
    m = Matching(tree1,tree2,node_cost)
    val = m.match()
    print('match =',val)
    assert val == 118
    print('distance table')
    print(m.getDistanceTable())


test_matching()