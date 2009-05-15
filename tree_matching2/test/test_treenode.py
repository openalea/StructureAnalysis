from openalea.tree_matching import *


def test_treenode():
    t = TreeNode(1,-1)
    assert t.father == -1,  "TreeNode member initialisation failed"

def test_treegraph():
    tree = TreeGraph()
    tree.addNode(0,-1)
    tree.addNode(1,0)
    tree.addNode(2,1)
    tree.addNode(3,1)
    tree.addNode(4,0)
    print "init"
    for i in range(5):
        print tree.getNode(i)    
    node = tree.getNode(4)
    node.addValue(2)
    print node
    print tree.getNode(4)


test_treenode()
test_treegraph()



def test_treenode_exists():
    tree = TreeGraph()
    tree.addNode(0,-1)
    n = tree.getNode(0)
    del tree
    n.addValue(1)
    print n

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

    print "#############"
    print tree1
    
    print "#############"
    print tree2

    node_cost = NodeCost()
    
    m = Matching(tree1,tree1,node_cost)
    m.match()
    print m.getDistanceTable()
    print m.getDBT(0,0)
    matrix = []
    for i in range(5):
        l = []
        for j in range(5):
            l += [m.getDBT(i,j)]
        matrix += [l]
    for i in range(5):
        print matrix[i]

         
    
#test_treenode_exists()
test_matching()

#from openalea.mtg.mtg import *
#import openalea.mtg.aml as wrap
#import openalea.aml as aml

#g = MTG("~/Documents/Travail/Aml/Pommier/B1001.mtg")

