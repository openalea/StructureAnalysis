from openalea.tree_matching import *


class MyNodeCost (NodeCost):
  def __init__(self):
     NodeCost.__init__(self)
  def getInsertionCost(self,a):
    return 18
  def getDeletionCost(self,a):
    return 50
  def getChangingCost(self,a,b):    
    #if a.id != b.id : return 1
    return 0
  def getMergingCost(self,l,b):
    return 3
    


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
  tree1 = TreeGraph()
  tree1.addNode(0)
  tree1.addNode(1,0)
  tree1.addNode(2,1)
  tree1.addNode(3,2)
  tree1.addNode(4,3)
  tree1.addNode(5,3)
  tree1.addNode(6,2)
  tree1.addNode(7,6)
  return tree1


def test_matching():
  tree1 = create_treegraph1()
  tree2 = create_treegraph2()

  node_cost = MyNodeCost()
  # Standard option = 0
  m = Matching(tree1,tree1,node_cost,0)
  m.verbose = False
  # Compact option = 1
  val = m.match()
  print('match =',val)
  print(m.getList(0,0))
  #for i in xrange(6):
  #  print m.getList(i,0)[0]
  #print
  #for i in xrange(6):
  #  print m.getList(0,i)[0]

def test_extmatching():
  tree1 = create_treegraph1()
  tree2 = create_treegraph2()

  node_cost = MyNodeCost()
  # Standard option = 0
  m = ExtMatching(tree1,tree2,node_cost,0)
  m.verbose = False
  # Compact option = 1
  val = m.match()
  print('match =',val)
  print(m.getList(0,0))
  #for i in xrange(6):
  #  print m.getList(i,0)[0]
  #print
  #for i in xrange(6):
  #  print m.getList(0,i)[0]


test_matching()
print("end Matching")
test_extmatching()

