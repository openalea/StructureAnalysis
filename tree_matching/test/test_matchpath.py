from openalea.tree_matching import *

class MyMatchPath (MatchPath):
  def __init__(self,l1,l2):
     MatchPath.__init__(self,l1,l2)
  def edgeCost(self,a,b):
#      print a,b,abs(a-b)
      return abs(a-b) 
  def minCostflow():
      return MatchPath.bipartiteMatching()
    

def test_createMatchpath():
    l1 = [10,20,30,60]
    l2 = [40,50]
    m = MyMatchPath(l1,l2)
    l = m.bipartiteMatching()
    print('**',l)



test_createMatchpath()

