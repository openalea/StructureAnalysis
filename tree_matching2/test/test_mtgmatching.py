from openalea.tree_matching.mtgmatching import *
from openalea.mtg import *

def my_random_mtg(nbvertices1 = 10,nbvertices2 = 10):
    g = MTG()
    # creating 2 scales
    root = g.root
    root1 = g.add_component(root)
    root1 = g.add_component(root1)
    
    # creating 
    vid = random_tree(g, root1, nb_vertices=nbvertices2)
    
    for i in xrange(nbvertices1-1):
        v1, complex1 = g.add_child_and_complex(vid)
        vid = random_tree(g, v1, nb_vertices=nbvertices2)

    return g

class MtgNodeCost(NodeCost):
    def __init__(self): NodeCost.__init__(self)
    def getDeletionCost(self,a) : return 1
    def getInsertionCost(self,b) : return 1
    def getChangingCost(self,a,b) : return 0
    
def test_mtgmatching():
    mtg1 = my_random_mtg()
    mtg2 = my_random_mtg()
    cost = MtgNodeCost()
    m = MtgMatching(mtg1,mtg2,scale1=2,scale2=2,cost=cost)
    m.match()
    matching = m.getList()
    totcost = m.getDBT()

class MtgExtNodeCost(NodeCost):
    def __init__(self): NodeCost.__init__(self)
    def getDeletionCost(self,a) : return 1
    def getInsertionCost(self,b) : return 1
    def getChangingCost(self,a,b) : return 0
    def getMergingCost(self,a,b) : return 0.1
    def getSplittingCost(self,a,b) : return 0.1

def test_mtgextmatching():
    mtg1 = my_random_mtg(5,5)
    mtg2 = my_random_mtg(5,5)
    cost = MtgExtNodeCost()
    m = MtgExtMatching(mtg1,mtg2,scale1=2,scale2=2,cost=cost)
    m.match()
    matching = m.getList()
    totcost = m.getDBT()

if __name__ == '__main__':
    test_mtgmatching()

