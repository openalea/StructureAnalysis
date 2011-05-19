from openalea.plantgl.all import *
from openalea.tree_matching import *
from vplants.pointreconstruction.util import *

class MyNodeCost (NodeCost):
  def __init__(self):
     NodeCost.__init__(self)
     
  def getInsertionCost(self,a):
	return 1
  
  def getDeletionCost(self,a):
    return 1
    
  def getChangingCost(self,a,b):    
	return 0
 
  def getMergingCost(self,a,b):
    return 2
    
  def getSplittingCost(self,a,b):
    return 0
	

def create_treegraph1():
  tree1 = TreeGraph()
  tree1.addNode(0)
  tree1.addNode(1,0)
  tree1.addNode(2,1)
  tree1.addNode(3,1)
  tree1.addNode(4,0)
  tree1.addNode(5,4)
  return tree1

def create_treegraph2():
  tree2 = TreeGraph()
  tree2.addNode(0)
  tree2.addNode(1,0)
  tree2.addNode(2,1)
  tree2.addNode(3,2)
  tree2.addNode(4,3)
  tree2.addNode(5,3)
  tree2.addNode(6,2)
  tree2.addNode(7,6)
  tree2.addNode(8,7)
  tree2.addNode(9,8)
  tree2.addNode(10,9)
  return tree2


def test_extmatching():
  tree1 = create_treegraph1()
  nbVertex1 = tree1.getNbVertex()
  tree2 = create_treegraph2()
  nbVertex2 = tree2.getNbVertex()

  node_cost = MyNodeCost()
  m = ExtMatching(tree1,tree2,node_cost,0)
  val = m.match()

  resultmatchinglist = m.getList(0,0)
  
  for m1,m2,mv in resultmatchinglist:
    print m1,m2,mv
    assert m1 < nbVertex1
    assert m2 < nbVertex2

if __name__ == '__main__':
    test_extmatching()
