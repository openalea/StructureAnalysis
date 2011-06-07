from openalea.tree_matching import *

class MyMatchPath (MatchPath):
  def __init__(self,l1,l2):
     MatchPath.__init__(self,l1,l2)
  def edgeCost(self,a,b):
      print a,b,abs(a-b)
      return abs(a-b) 
    

def test_createMatchpath():
    l1 = [1,2,3]
    l2 = [4,5]
    m = MyMatchPath(l1,l2)
    l = m.minCostFlow()
    print '**',l



test_createMatchpath()

