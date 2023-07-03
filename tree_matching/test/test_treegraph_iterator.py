from openalea.tree_matching import *

def return_self(self): return self


def create_treegraph():
	treegraph = TreeGraph()
	treegraph.addNode(0)
	treegraph.addNode(1,0)
	treegraph.addNode(2,1)
	treegraph.addNode(3,2)
	treegraph.addNode(4,3)
	treegraph.addNode(5,3)
	treegraph.addNode(6,1)
	treegraph.addNode(7,6)
	treegraph.addNode(8,7)
	treegraph.addNode(9,8)
	treegraph.addNode(10,9)
	return treegraph

t = create_treegraph()

def iterator_check(r = 2):
	l1 = list(t.subtree_iter(r))
	l2 = t.subtree_list(r)
	print(l1)
	assert l1 == l2

def test_iterator():
   for i in range(10):
       yield iterator_check, i
	
if __name__ == '__main__':
	test_iterator()