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


def test_choicetable():
  tree1 = create_treegraph1()
  tree2 = create_treegraph2()

  node_cost = MyNodeCost()
  m = Matching(tree1,tree2,node_cost,0)
  val = m.match()
  a = m.getList(0,0)
  a = [(i,j) for i,j,k in a]

  # Get the choice table from matching  
  ct = m.getChoiceTable()
  # test for similar getList matching value
  b = ct.getList(0,0,tree1,tree2)
  assert a == b

  # dump and load a choice table  
  fname = 'test_choicetable.dump'
  ct.dump(fname)
  ct2 = ChoiceTable.load(fname)
  
  # check for getList to work.
  b = ct2.getList(0,0,tree1,tree2)
  assert a == b
  
if __name__ == '__main__':
    test_choicetable()

