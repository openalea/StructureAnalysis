# a test for the class etrees.Tree: constructor and basic methods
import sys, os
import openalea.tree_statistic.trees.etrees as etrees
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tool

inf_bound=0
sup_bound=3
probability= 0.6
ident=stat_tool.DistributionIdentifier.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool._ParametricModel(ident, inf_bound, sup_bound, 
                                    parameter, probability)
print distrib
max_depth= 3
max_size= 8
n=5
T=trees.etrees.Tree([0], 0, 0)
v=range(n)
for i in range(n):
    v[i]=T.AddVertex([i+1])
T.AddEdge(v[0], v[1])
T.AddEdge(v[0], v[2], '<')
T.AddEdge(v[1], v[3])
T.AddEdge(v[1], v[4], '<')
print "The root of T is ", T.Root()
print "The parent of vertex 4 is:", T.Parent(4)
print "A predefined tree: "
print T
print "... with depth", T.Depth(), " and size ", T.Size(), "."
for i in range(n):
    print "Number of children of node ", i+1, ": ", T.NbChildren(i)
# Extract a subtree
print "Subtree of T rooted at node 2:" 
print T.Display(False, True, 1)
# tree iterators
vertex_iterator=iter(T)
vertices=[]
try:
    while 1:
        vertices.append(vertex_iterator.next()+1)
except StopIteration: pass
print "A vertex iterator:", vertices
preorder=T.Preorder()
vertices=[]
try:
    while 1:
        vertices.append(preorder.next()+1)
except StopIteration: pass
print "A preorder tree traversing:", vertices
inorder=T.Inorder()
vertices=[]
try:
    while 1:
        vertices.append(inorder.next()+1)
except StopIteration: pass
print "An inorder tree traversing:", vertices
postorder=T.Postorder()
vertices=[]
try:
    while 1:
        vertices.append(postorder.next()+1)
except StopIteration: pass
print "An postorder tree traversing:", vertices
breadthorder=T.Breadthorder()
vertices=[]
try:
    while 1:
        vertices.append(breadthorder.next()+1)
except StopIteration: pass
print "A breadth-fisrt tree traversing:", vertices
leavesfirstorder=T.LeavesFirstorder()
vertices=[]
try:
    while 1:
        vertices.append(leavesfirstorder.next()+1)
except StopIteration: pass
print "A tree traversing that minimizes the distance " \
      "to nearest leaf node:", vertices
# select subtree
ST=T.SelectSubTree(1)
print "Extract subtree rooted at vertex 2 from T:"
print ST
ST=T.SelectSubTree(1, False)
print "Prune subtree rooted at vertex 2 from T:"
print ST
# extract sequences
import openalea.aml as amlPy
TS=trees.Trees([T])
S1=TS.BuildSequences(True)
S2=TS.BuildSequences(False)
amlPy.Display(S1, ViewPoint="Data")
amlPy.Display(S2, ViewPoint="Data")
# W=T creates an alias
W=trees.etrees.Tree(T)
W.Put(W.Root(), [100])
print "A (non-alias) copy of T with modified root attribute:"
print W
##C=trees.etrees.TreeStructure(0, 1)
##print "A default tree structure", C
##v=range(n)
##for i in range(n):
##    v[i]=C.AddVertex()
##C.AddEdge(v[0], v[1])
##C.AddEdge(v[0], v[2], '<')
##C.AddEdge(v[1], v[3])
##C.AddEdge(v[1], v[4], '<')
##print "A tree structure build manually:"
##C.Display()
##print "... with depth", C.Depth(), " and size ", C.Size(), "."
##C=trees.etrees.TreeStructure(T)
##print "Extracting a structure from T:"
##print C
print "A random tree structure R:"
R=trees.etrees.TreeStructure(distrib, max_size, max_depth)
print R
print "(size, depth) = (", R.Size(), ", ", R.Depth(), ")."
print "Copying this structure to T:"
T=trees.etrees.Tree([1.0, 0], R)
print T
print "Random simulation of the attributes of T:"
distrib_list=[]
for i in range(T.NbInt()):
    distrib_list.append(distrib)
T.Simulate(distrib_list)
print T
print "Conversion from etrees to tree"
PT=trees.Tree(T)
print PT
T=trees.etrees.Tree(PT)
print T
# S()
try:
    print "The parent of the root vertex is:", T.Parent(T.Root())
except IndexError, e:
    print e
else:
    print "Failed to detect incorrect call to Parent"

#commented those lines that raise error message. Author(s) should fix it. TC.Dec 2008. 
#if (T.IsRoot(4)):
#    raise WarningError, "Failed to detect that vertex 4 is not the root"
#else:
#    print "Is vertex 4 the root ? ", str(T.IsRoot(4))
