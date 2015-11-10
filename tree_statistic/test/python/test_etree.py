# -*- coding: utf-8 -*-
# a test for the class etrees.Tree: constructor and basic methods
import sys, os
import openalea.tree_statistic.trees.etrees as etrees
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tool

def init():
    """Defines some constants"""
    global distrib, max_depth, max_size, nbtrees, n
    stat_tool.plot.DISABLE_PLOT = True
    inf_bound = 1
    sup_bound = 3
    distrib = stat_tool.Uniform(inf_bound, sup_bound)
    nbtrees = 40
    max_depth = 3
    max_size = 8
    n = 5

init()

# Build a predefined tree
T = trees.etrees.Tree([0], 0, 0)
v = range(n)
for i in range(n):
    v[i] = T.AddVertex([i+1])
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
vertex_iterator = iter(T)
vertices = []
try:
    while 1:
        vertices.append(vertex_iterator.next()+1)
except StopIteration: pass
print "A vertex iterator:", vertices
preorder = T.Preorder()
vertices = []
try:
    while 1:
        vertices.append(preorder.next()+1)
except StopIteration: pass
print "A preorder tree traversing:", vertices
inorder = T.Inorder()
vertices = []
try:
    while 1:
        vertices.append(inorder.next()+1)
except StopIteration: pass
print "An inorder tree traversing:", vertices
postorder = T.Postorder()
vertices = []
try:
    while 1:
        vertices.append(postorder.next()+1)
except StopIteration: pass
print "An postorder tree traversing:", vertices
breadthorder = T.Breadthorder()
vertices = []
try:
    while 1:
        vertices.append(breadthorder.next()+1)
except StopIteration: pass
print "A breadth-fisrt tree traversing:", vertices
leavesfirstorder = T.LeavesFirstorder()
vertices = []
try:
    while 1:
        vertices.append(leavesfirstorder.next()+1)
except StopIteration: pass
print "A tree traversing that minimizes the distance " \
      "to nearest leaf node:", vertices
# select subtree
ST = T.SelectSubTree(1)
print "Extract subtree rooted at vertex 2 from T:"
print ST
ST = T.SelectSubTree(1, False)
print "Prune subtree rooted at vertex 2 from T:"
print ST
# extract sequences
import openalea.aml as amlPy
TS = trees.Trees([T])
S1 = TS.BuildSequences(True)
S2 = TS.BuildSequences(False)
amlPy.Display(S1, ViewPoint="Data")
amlPy.Display(S2, ViewPoint="Data")
# W=T creates an alias
W = trees.etrees.Tree(T)
W.Put(W.Root(), [100])
print "A (non-alias) copy of T with modified root attribute:"
print W

# Tree structures
print "A random tree structure R:"
R = trees.etrees.TreeStructure(distrib, max_size, max_depth)
print R
print "(size, depth) = (", R.Size(), ", ", R.Depth(), ")."
# Copy constructor
R2 = trees.etrees.TreeStructure(R)
# Add vertex and egde
v = R2.AddVertex()
# should work but does not ?
#R2.AddEdge(2, v)
#print R2.Display()

print "Copying this structure to T:"
T = trees.etrees.Tree([1.0, 0], R)
print T
print "Random simulation of the attributes of T:"
distrib_list = []
for i in range(T.NbInt()):
    distrib_list.append(distrib)
T.Simulate(distrib_list)
print T
print "Conversion from etrees to tree"
PT = trees.Tree(T)
print PT
T = trees.etrees.Tree(PT)
print T
# S()
try:
    print "The parent of the root vertex is:", T.Parent(T.Root())
except IndexError, e:
    print e
else:
    print "Failed to detect incorrect call to Parent"

if (T.IsRoot(4)):
   raise WarningError, "Failed to detect that vertex 4 is not the root"
else:
   print "Is vertex 4 the root ? ", str(T.IsRoot(4))

# test copy and setting tree vids
print "Read Trees object from MTG file 'sample_mtg_forest.txt':"
T = trees.etrees.Trees("sample_mtg_forest.txt")
print "Read ", T.NbTrees(), " trees."
TreeVidDictList = []
TreeVidDictListCp = []
for t in range(T.NbTrees()):
    TreeVidDictList.append(T.TreeVertexId(t))
    TreeVidDictListCp.append(T.TreeVertexId(t))
    for k in TreeVidDictListCp[t].keys():
        TreeVidDictListCp[t][k] += 1

T._SetMTGVidDictionary(TreeVidDictListCp)
TreeVidDictListCp = []
for t in range(T.NbTrees()):
    TreeVidDictListCp.append(T.TreeVertexId(t))

try:
    T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
except IndexError, msg:
    print msg
else:
    assert False, "Failed to raise IndexError"

TreeVidDictListCp = []
for t in range(T.NbTrees()):
    TreeVidDictListCp.append(dict(TreeVidDictList[t]))

TreeVidDictListCp[0][3] = 0

try:
    T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
except ValueError, msg:
    print msg
else:
    assert False, "Failed to raise ValueError"

TreeVidDictListCp = []
for t in range(T.NbTrees()):
    TreeVidDictListCp.append(dict(TreeVidDictList[t]))
    for k in TreeVidDictListCp[t].keys():
        TreeVidDictListCp[t][k] = str(TreeVidDictListCp[t][k])


T._SetMTGVidDictionary(TreeVidDictListCp)
try :
    T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
except TypeError, msg:
    print msg
else:
    assert False, "Failed to raise TypeError"


TreeVidDictListCp = []
for t in range(T.NbTrees()):
    TreeVidDictListCp.append({})
    for k in TreeVidDictList[t].keys():
        TreeVidDictListCp[t][str(k)] = TreeVidDictList[t][k]

T._SetMTGVidDictionary(TreeVidDictListCp)
try :
    T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
except Exception, msg:
    print msg
else:
    assert False, "Failed to raise TypeError-like Exception"


