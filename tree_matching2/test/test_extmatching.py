from openalea.tree_matching import *


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
  nbVertex1 = tree1.getNbVertex() # 6
  tree2 = create_treegraph2()
  nbVertex2 = tree2.getNbVertex() # 11

  node_cost = MyNodeCost()
  m = ExtMatching(tree1,tree2,node_cost,0)
  val = m.match()

  resultmatchinglist = m.getList(0,0)
  
  # we suppose that because split is 0 cost and merge is big cost, no merging will occur and element of tree2 will be match only one time
  m2list = set() 
  for m1,m2,mv in resultmatchinglist:
    print m1,m2,mv
    # check if it is a valid element of tree 1
    assert 0 <= m1 < nbVertex1
    # check if it is a valid element of tree 2
    assert 0 <= m2 < nbVertex2
    # check if there is no merging
    assert not m2 in m2list
    m2list.add(m2)
    
"""
An error occur that found an extra node with strange id at the last line for m2
Here comes the log
**********
0 2 -1.0
0 1 -1.0
0 0 -1.0
1 3 0.0
2 5 0.0
3 4 0.0
4 6 0.0
5 9 -1.0
5 8 -1.0
5 7 -1.0
5 52443312 -1.0
--------------
AssertionError
...
**********
On windows the strange id is usually a number between -1 and 10
"""

if __name__ == '__main__':
    test_extmatching()
