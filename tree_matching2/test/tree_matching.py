from openalea.plantgl.all import *
from openalea.tree_matching import *


class MyNodeCost (NodeCost):
  def __init__(self):
     NodeCost.__init__(self)
  def getInsertionCost(self,a):
    return 18
  
  def getDeletionCost(self,a):
    return 50
    
  def getChangingCost(self,a,b):    
    if a.id != b.id : return 1
    return 0


def create_treegraph():
  tree1 = TreeGraph()
  tree1.addNode(0)
  tree1.addNode(1,0)
  tree1.addNode(2,1)
  tree1.addNode(3,1)
  tree1.addNode(4,0)
  tree1.addNode(5,4)
  return tree1


def test_matching():
  tree1 = create_treegraph()
  tree2 = create_treegraph()

  node_cost = MyNodeCost()
  m = Matching(tree1,tree2,node_cost)
  val = m.match()
  print 'match =',val
  for i in xrange(6):
    print m.getList(i,0)[0]
  print
  for i in xrange(6):
    print m.getList(0,i)[0]


test_matching()

