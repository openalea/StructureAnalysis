import sys
import os
import string
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.int_fl_containers as int_fl_containers
#from openalea.tree_statistic.trees import *

import openalea.tree_statistic.trees.etrees as etrees

ctree = etrees.ctree



# TreeValues
print "Building a TreeValue [2, 1.0, 0]"
default_value = int_fl_containers.Int_fl_container(1, 0)
i = ctree.TreeValue([0, 0., 0])
i[0] = 1
i[1] = 1
i[0] = 2
print "i: ", i
ifc = int_fl_containers.Int_fl_container(i)
print "Int_fl_container conversion of i: ",ifc


# test of a default tree
print "Test of a default CTree:"
w = ctree.CTree(1,0,0,1)

# nb int, nb float, root, nb vertices
print "Is vertex 0 the root vertex ?: ", w.IsRoot(0)
print "Is vertex 1 the root vertex ?: ", w.IsRoot(1)
print "Root vertex: ", w.Root()
print "Tree size: ", w.Size()
print "Entire tree: ", w
    
# test of a manually-defined tree
print "Test of a manually-defined tree"
tm=ctree.CTree(1,0,0,0)
vlist=[]
tm_size=5
for i in range(5):
    value=ctree.TreeValue([i])
    vlist.append(tm.AddVertex(value))
tm.AddEdge(vlist[0], vlist[1], True)
tm.AddEdge(vlist[0], vlist[2], False)
tm.AddEdge(vlist[1], vlist[3], True)
tm.AddEdge(vlist[1], vlist[4], True)
print tm
for i in range(5):
    if not(tm.IsRoot(i)):
        print "Parent of vertex ", i, ": ", tm.Parent(i)
        print "Is edge (", tm.Parent(i), ", ",i,") a valid edge ?:"
        print tm.IsEdge(tm.Parent(i), i)

# test of the copy constructor
print 'Test of copy constructor: a copy of the above tree'
cptm = ctree.CTree(tm)
print cptm
for i in range(5):
    if not(tm.IsRoot(i)):
        print "Parent of vertex ", i, ": ", cptm.Parent(i)
        print "Is edge (", cptm.Parent(i), ", ",i,") a valid edge ?:", \
        cptm.IsEdge(cptm.Parent(i), i)
inf_bound=0
sup_bound=3
probability= 0.6
ident=stat_tool.DistributionIdentifierType.UNIFORM
parameter=stat_tool.D_DEFAULT
distrib= stat_tool.distribution._DiscreteParametricModel(ident, inf_bound, sup_bound, parameter, probability)

print distrib
max_depth= 3
max_size= 20
print "Default TreeStructure object"
u=ctree.TreeStructure()
print u
print "Defining a TreeStructure object"
u=ctree.TreeStructure(0, 1)

#with maximal depth ", max_depth, "."

u.Simulate(distrib, max_size, max_depth)
print u # u.Display(u.Root())
# is_edge
print "Is (0, 1) an edge ? ", u.IsEdge(0, 1)
if (u.IsEdge(0, 1)):
    print "Type of this edge: ", u.EdgeType(0, 1)
print "Is (0, ",sup_bound+1,") an edge ? ", u.IsEdge(0, sup_bound+1)
if (u.IsEdge(0, sup_bound+1)):
    print "Type of this edge: ", u.EdgeType(0, sup_bound+1)
print "Effective size: ", u.Size()
print "Creating directly a new tree with random structure"
u2= ctree.TreeStructure(distrib, max_size, max_depth)
print u2 # u2.Display(u2.Root())
# select subtree
print "Size of u2:", u2.Size()
ST=ctree.SelectSubTreeStructure(u2, 3, True)
print "Extracting subtree rooted at vertex 3 from above tree:"
print ST
ST=ctree.SelectSubTreeStructure(u2, 3, False)
print "Pruning subtree rooted at vertex 3 from above tree:"
print ST
print "Building a tree with attribute from random structure u2"
print "(attributes are vertex ids): "
w= ctree.CTree(u2, default_value)
tree_iterator= iter(w)
w_vertices=[]
try:
    while 1:
        w_vertices.append(tree_iterator.next())
except StopIteration: pass
for i in w_vertices:
    value=ctree.TreeValue([i])
    w.Put(i, value)
print w
for i in range(w.Size()):
    if not(w.IsRoot(i)):
        print "Parent of vertex ", i, ": ", w.Parent(i)
        print "Is edge (", w.Parent(i), ", ",i,") a valid edge ?:", \
        w.IsEdge(w.Parent(i), i)
u=w.GetStructure()
print "Copying structure to w and u:"
print u 
u.Display()
# select subtree
print "Extracting subtree rooted at vertex 3 from w:"
ST=ctree.SelectSubTree(w, 3, True)
print ST
print "Pruning subtree rooted at vertex 3 from w:"
ST=ctree.SelectSubTree(w, 3, False)
print ST
emptyv=ctree.TreeValue([w.Size()])
# AddVertex
print "Add vertex ", w.Size()
v=w.AddVertex(emptyv)
w.AddEdge(w.Root(), v, False)
print w # w.Display(w.Root())
# Copy a structure
print "Simulate and copy a tree structure: "
u.Simulate(distrib, max_size, max_depth)
i= w.Get(w.Root())
print "Value of the root vertex: ", i
w.Put(w.Root(), emptyv)
w.SetStructure(u, i)
print w
print "Number of integer variables: ", w.NbInt()
# Tree iterators
print "Building a vertex iterator for w: "
tree_iterator= w.__iter__()
vertices=[]
try:
    while 1:
        vertices.append(tree_iterator.next())
except StopIteration: pass
print vertices
print "Check the size compatibility. w:", w.Size(), " / iterator: ", len(vertices)
tree_iterator= w.Vertices()
vertices=[]
# Children iterators
print "Building a children iterator for w at root vertex:"
children_iterator= w.Children(w.Root())
children=[]
try:
    while 1:
        children.append(children_iterator.next())
except StopIteration: pass
print "Check the compatibility for number of children.", \
"NbChildren: ", w.NbChildren(w.Root()), " / iterator: ", len(children)
print "Children of the root of w: ", children
tree_iterator= u2.__iter__()
vertices=[]
try:
    while 1:
        vertices.append(tree_iterator.next())
except StopIteration: pass
print "Size of u2:", u2.Size(), " == ", len(vertices)
print "Vertices of u2: ", vertices
# Building a ctree from a MTG
print "Building a tree from a MTG"
from amlPy import *
if not getmode():
    s=1
    setmode(1)
# conversion from AML object to Python
def f1(x):
    return Feature(x,"length")
def f2(x):
    return Feature(x,"fruits")
M=MTG("sample_mtg.txt")
nbfloat=0
nbint=0
mtg_file=file("sample_mtg.txt", 'r')
mtg_code=mtg_file.readline()
while (mtg_code.upper().find("FEATURES")==-1):
    mtg_code=mtg_file.readline()
while (string.find(string.upper(mtg_code),"NAME")==-1):
    mtg_code=mtg_file.readline()
mtg_code=mtg_file.readline()
while mtg_code.isspace():
    mtg_code=mtg_file.readline()
s=0
r=MTGRoot()
Tmtg=ctree.CTree(1,1,r,0)
_v=ComponentRoots(r, Scale=1)
print "_v: ", _v
_v=_v[0]
mtg2tree_vertices={}
_b=Descendants(_v) # warning : contains _v !
print "Size", Tmtg.Size()
attribute_list=[f1,f2]
for _x in _b:
    values=[f1(_x), f2(_x)]
    mtg2tree_vertices[_x]=Tmtg.AddVertex(ctree.TreeValue(values))
    if (_x != _v):
        Tmtg.AddEdge(mtg2tree_vertices[Father(_x)],mtg2tree_vertices[_x], False)
print Tmtg
if s:
    setmode(0)
