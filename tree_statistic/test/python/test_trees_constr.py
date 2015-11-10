# -*- coding: utf-8 -*-
# a test for the class trees.Trees: constructor and basic methods
import sys
import os
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tools


from PyQt4.QtCore import *
from PyQt4.QtGui import *


class MyThread(QThread):
    def __init__(self):
        QThread.__init__(self)
    
    def doit(self, T):
        app = QApplication([])
        T.Plot()
        #self.start()
        #app.exec_()


stat_tools.plot.DISABLE_PLOT = True

inf_bound=0
sup_bound=3
distrib = stat_tools.Uniform(inf_bound, sup_bound)
print "Distribution used for the number of children and the tree attributes:"
print distrib
max_depth=3
max_size=10
nbtrees=40
# define a set of trees
tree_list=[]
tv=[1., 0, 1, 2.] # trees.TreeValue([1., 0])
R=trees.TreeStructure(distrib, max_size, max_depth)
tmp_tree=trees.Tree(tv, R)
n=1
tree_list.append(trees.Tree(tmp_tree))
while n < nbtrees:
    n=n+1
    R.Simulate(distrib, max_size, max_depth)
    tmp_tree=trees.Tree(tv, R)
    tree_list.append(trees.Tree(tmp_tree))    
distrib_list=[]
for i in range(tmp_tree.NbInt()):
    distrib_list.append(distrib)
    
for n in range(len(tree_list)):
    tree_list[n].Simulate(distrib_list)
# initialize a Trees object
T=trees.Trees(tree_list)
print "A set of trees with variable types", T.Types(), \
      " (initialized by ", tv, "):"
print T
print "Display Trees object with detail level 2:"
T.Display(Detail=2)
print "Display Trees object with data and detail level 1:"

T.Display(ViewPoint="DATA", Detail=1)
T.ExtractHistogram("Value", variable=1).plot()
T.Plot(ViewPoint="FirstOccurrenceRoot", variable=1)
T.MPlot(ViewPoint="FirstOccurrenceRoot", variable=1)

# save trees as openalea.mtg.MTG
T.PickleDump("tree.pkl")
T2 = trees.PickleLoad("tree.pkl")
#T.Plot()

#m = MyThread()
#m.doit(T)

print "Number of trees (supposedly ", nbtrees, "): ", T.NbTrees()
print "Comparison of trees contained in Trees " \
      "and those used for initialization: "
equal=True
differences=[]
for n in range(len(tree_list)):
    if str(T.Tree(n))!=str(tree_list[n]):
        equal=False
        differences.append(n)
if not equal:
    print "Differences found: ", differences
else:
    print "Match!"
    # print "Trees.Tree(", n, "):"
    # print T.Tree(n)
    # print "Tree number ", n, ":"
    # print tree_list[n]
    
# trees with one vertex
print "Test using trees of size one"
omax_size=1
onbtrees=2
# define a set of trees
otree_list=[]
otv=[1., 0, 1, 2.] # trees.TreeValue([1., 0])
OR=trees.TreeStructure(distrib, omax_size, max_depth)
otmp_tree=trees.Tree(otv, OR)
on=1
otree_list.append(trees.Tree(otmp_tree))

while on < onbtrees:
    on=on+1
    OR.Simulate(distrib, omax_size, max_depth)
    otmp_tree=trees.Tree(otv, OR)
    otree_list.append(trees.Tree(otmp_tree))    
for on in range(len(otree_list)):
    otree_list[on].Simulate(distrib_list)
# initialize a Trees object
OT=trees.Trees(otree_list)
print "A set of trees with variable types", T.Types(), \
      " (initialized by ", tv, "):"
print OT
# copy constructor
print "Test of the constructor by copy: "
T2=trees.Trees(T)
print T2
# constructor from a MTG
print "Read Trees object from MTG file 'sample_mtg_forest.txt':"
T=trees.Trees("sample_mtg_forest.txt")
print "Read ", T.NbTrees(), " trees."
print T
print "Print 1st tree of Trees object:"
T.Tree(0).Display()



T=T.Cluster("Step", 0, 10)


import openalea.aml as amlPy
amlPy.setmode(1)
# build vectors from trees
vec = T.BuildVectors()
print vec
# reading a Tree from a MTG with filter and custom attributes
print "Read tree from MTG file 'sample_mtg_forest.txt' using a filter", \
      "and custom attributes:"
f=lambda x: amlPy.Feature(x, "Length")*3*(amlPy.Feature(x, "Diam")/2)**2
# FM=amlPy.MTG("sample_mtg_forest.txt")
# f=lambda x: amlPy.Feature(x, "Diam")
# print [amlPy.Feature(2, "Length"), amlPy.Feature(2, "Diam")]
# print f(5)
filter=lambda x: x < 6
# filter=lambda x: True
attributes=["anything"]
T=trees.Trees("sample_mtg_forest.txt", filter, attributes, [f], scale=2)
print "Read ", T.NbTrees(), " trees."
print T
print "Print 1st tree of Trees object:"
print T.Tree(0).Display()
# checking exceptions
# invalid file name


try:
    T=trees.Trees("no_such_file.t\/t")
except IOError, e:
    print e
else:
    print "Failed to raise exception for invalid MTG file name"
    raise e
# inconsistent attribute number
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [], [f], scale=2)
except ValueError, v:
    print v
else:
    print "Failed to raise exception for inconsistent attribute number"
# bad attribute name

try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [f], [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad attribute function
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes, attributes,
                  scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# bad attribute type
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes,
    [lambda x: amlPy.Descendants(x)], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
# bad attribute name
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, [f], [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad filtering function
try:
    T=trees.Trees("sample_mtg_forest.txt", attributes, attributes, [f], scale=2)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# filtering function does not filter descendants
try:
    T=trees.Trees("sample_mtg_forest.txt", lambda x: x != 2, attributes, [f],
                  scale=2)
except IndexError, i:
    print i
else:
    print "Failed to raise exception for filter not filtering descendants"
try:
    T=trees.Trees("sample_mtg_forest.txt", filter, attributes, [f], scale=1)
except TypeError, t:
    print t

# build MTG
import openalea.mtg as mtg
vtxlist = [2*i for i in range(2,27)]
mtgvtx = []
g = mtg.MTG()
mtgvtx += [g.add_component(0)]
mtgvtx += [g.add_component(0)]
mtgroots = list(mtgvtx)
for v in range(13):
    mtgvtx += [g.add_component(mtgroots[0], vtxlist[v])]

g.add_child(4,6,edge_type='<')
for v in range(2):
    g.add_child(4,8+2*v,edge_type='+')
g.add_child(6,12,edge_type='<')
for v in range(2):
    g.add_child(6,14+2*v,edge_type='+')
for v in range(4):
    g.add_child(8,18+2*v,edge_type='+')
g.add_child(24,26,edge_type='<')
g.add_child(24,28,edge_type='+')

for v in range(13, 25):
    mtgvtx += [g.add_component(mtgroots[1], vtxlist[v])]
g.add_child(30,32,edge_type='<')
for v in range(3):
    g.add_child(30,34+2*v,edge_type='+')
g.add_child(34,40,edge_type='<')
for v in range(2):
    g.add_child(34,42+2*v,edge_type='+')
for v in range(4):
    g.add_child(44,46+2*v,edge_type='+')
g.add_property("Name1")
g.add_property("Name2")
for v in mtgvtx:
    g.node(v).Name1 = v
    g.node(v).Name2 = 52-v
from openalea.mtg import treestats
flist = [lambda x: g.node(x).Name1, lambda x: g.node(x).Name2]
T = treestats.extract_trees(g, 2, lambda x: True, flist)
T.PickleDump("tree.pkl")
T2 = trees.PickleLoad("tree.pkl")
print T2.NbTrees()
T2.Tree(0).Display()

