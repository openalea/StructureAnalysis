from openalea.tree_matching import *

def create_treegraph():
    tree1 = TreeGraph()
    tree1.addNode(0,-1)
    tree1.addNode(1,0)
    tree1.addNode(2,1)
    tree1.addNode(3,1)
    tree1.addNode(4,0)
    return tree1


def test_treenode_persistency():
    tree = create_treegraph()
    n = tree.getNode(0)
    n.toto = True
    m = tree.getNode(0)
    if id(n) != id(m):
        import warnings
        warnings.warn('object identity not preserved in this case')
        

def test_treenode_persistency2():
    tree = TreeGraph()
    l = TreeNode(0)
    tree.addNode(l)
    n = tree.getNode(0)
    n.toto = True
    m = tree.getNode(0)
    assert m.toto
    assert id(l) == id(n) == id(m) and 'object identity not preserved even when creating TreeNode from python'
    
if __name__ == '__main__':
    test_treenode_persistency()
    test_treenode_persistency2()    