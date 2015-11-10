# -*- coding: utf-8 -*-
# a test for the class trees.trees.Trees: constructor and basic methods
 
# The classes in PyQt4 are useful to test
# visualization by PlantFrame

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import sys, os
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.trees as trees


class MyThread(QThread):
    def __init__(self):
        QThread.__init__(self)

    def doit(self, P):
        app = QApplication([])
        amlPy.Plot(P)
        #self.start()
        #app.exec_()


def init():
    """Defines some constants"""
    global distrib, max_depth, max_size, nbtrees, n
    stat_tool.plot.DISABLE_PLOT = True
    inf_bound = 1
    sup_bound = 3
    distrib = stat_tool.Uniform(inf_bound, sup_bound)
    nbtrees = 40
    max_depth = 3
    max_size = 20
    n = 5

init()
print distrib

################################
# TreeStructure
################################

D = trees.TreeStructure()
print "Default tree structure: ", D
print "A random tree structure R:"
R = trees.TreeStructure(distrib)
print R
print "(size, depth) = (", R.Size(), ", ", R.Depth(), ")."
R.Simulate(distrib, max_size, max_depth)
print "A second realization of R: "
print R
# TreeStructure.SelectSubTree
ST = R.SelectSubTree(1)
print "Extract subtree rooted at vertex 1 from R:"
print ST
# pruning
ST = R.SelectSubTree(1, False)
print "Prune subtree rooted at vertex 1 from R:"
print ST

# checking exceptions
try:
    ST = R.SelectSubTree(R.Size()+1)
except IndexError, i:
    print i
else:
    print "Failed to raise exception for bad vertex identifier"
try:
    ST = R.SelectSubTree(R)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad vertex identifier type"

################################
# Tree
################################
# TreeStructure to Tree
print "Copy the entire structure to a tree T:"
T = trees.Tree([1.0, 0], R)

T.Display(False)
print "Root of T:", T.Root()
print "Size of T:", T.Size()
print "Display T with vids: "
T.Display(True)
# Tree.SelectSubtree
ST = T.SelectSubTree(3)
print "Extract subtree rooted at vertex 3 from T:"
print ST
ST = T.SelectSubTree(3, False)
print "Prune subtree rooted at vertex 3 from T:"
print ST
# checking exceptions
try:
    ST = T.SelectSubTree(-1)
except IndexError, i:
    print i
else:
    print "Failed to raise exception for bad vertex identifier"
try:
    ST = T.SelectSubTree(T)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad vertex identifier type"
# Tree to TreeStructure
print "Extracting the structure from T:"
R = trees.TreeStructure(T)
print R
print "Display the vertex identifiers of the structure:"
R.Display()
# Tree.Simulate
print "Random simulation of the attributes of T and modification of its root:"
distrib_list = []
for i in range(T.NbInt()):
    distrib_list.append(distrib)
T.Simulate(distrib_list)
# Tree.EdgeType
print "Edge (0,1) has type: ", T.EdgeType(0,1)
# Tree.Put
T.Put(T.Root(), [1.8, 1])
T.Display(True)
# checking exceptions
try:
    T.Put(T.Root(), 1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
try:
    T.Put(T.Root(), "root")
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
try:
    T.Put(T.Root(), [1., 2.])
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
try:
    T.EdgeType(T.Size(), 0)
except IndexError, t:
    print t
else:
    print "Failed to raise exception for bad edge identifier"
try:
    T.EdgeType("root", 0)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad edge type"
# Tree.Round
print "Round values of a tree"
val = T.Get(2)
val[0] = val[0]+0.1
val = T.Put(2, val)
RT = T.Round()
print RT
# Tree.Save
print "Save tree T into file 'fake_mtg.txt'"
T.Save("fake_mtg.txt", True, ["longueur", "fruits"])
# checking exceptions
try:
    T.Save("fake_mtg.txt")
except IOError, i:
    print i
else:
    print "Failed ro raise exception for overwriting file"
try:
    T.Save("fake_mtg.txt", True, ["longueur"])
except IndexError, t:
    print t
else:
    print "Failed ro raise exception for bad number of attributes"
try:
    T.Save("fake_mtg.txt", True, [1, 2])
except TypeError, t:
    print t
else:
    print "Failed ro raise exception for attribute name type"
try:
    T.Save("fake_mtg.txt", True, ["var1", "bad name"])
except ValueError, v:
    print v
else:
    print "Failed ro raise exception for bad attribute name"



# reading a Tree from a MTG
print "Read tree from MTG file 'sample_mtg.txt':"
import openalea.aml as amlPy
amlPy.setmode(1)
MT = trees.Tree("sample_mtg.txt")
print MT.Display()
MT.Save("sample_mtg_save.txt", True)
# reading a Tree from a MTG with filter and custom attributes
print "Read tree from MTG file 'sample_mtg.txt' using a filter ", \
      "and custom attributes:"
f = lambda x: amlPy.Feature(x, "length")*amlPy.Feature(x, "fruits")
# filter=lambda x: x != 2
filter = lambda x: x < 6
attributes = ["anything"]
MT = trees.Tree("sample_mtg.txt", filter, attributes, [f], scale=1)
print MT.Display()


import openalea.aml as amlPy

MSTG = amlPy.MTG("sample_mtg.txt")
DR = amlPy.DressingData("dressing.drf")

def Longueur(x):
    return 20*amlPy.Feature(x, "length")
def Rayon(x):
    return 1

P = amlPy.PlantFrame(amlPy.MTGRoot(), Scale=1, DressingData = DR, Length = Longueur, 
                     BottomDiameter=Rayon)


# this is a call to the QThread, that is required to prevent QTapplication seg fault
m = MyThread()
m.doit(P)


# checking exceptions
# invalid file name
try:
    MT = trees.Tree("no_ $uch file.t\/t")
except IOError, e:
    print e
else:
    print "Failed to raise exception for invalid MTG file name"
    raise e
# inconsistent attribute number
try:
    MT = trees.Tree("sample_mtg.txt", filter, [], [f], scale=1)
except ValueError, v:
    print v
else:
    print "Failed to raise exception for inconsistent attribute number"
# bad attribute name
try:
    MT = trees.Tree("sample_mtg.txt", filter, [f], [f], scale=1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad attribute function
try:
    MT = trees.Tree("sample_mtg.txt", filter, attributes, attributes, scale=1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# bad attribute type
try:
    MT = trees.Tree("sample_mtg.txt", filter, attributes, 
    [lambda x: amlPy.Descendants(x)], scale=1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute type"
# bad attribute name
try:
    MT = trees.Tree("sample_mtg.txt", filter, [f], [f], scale=1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute name"
# bad filtering function
try:
    MT = trees.Tree("sample_mtg.txt", attributes, attributes, [f], scale=1)
except TypeError, t:
    print t
else:
    print "Failed to raise exception for bad attribute function"
# filtering function does not filter descendants
try:
    MT = trees.Tree("sample_mtg.txt", lambda x: x != 2, attributes, [f], scale=1)
except IndexError, i:
    print i
else:
    print "Failed to raise exception for filter not filtering descendants"
# reading a forest MTG
print "Read a forest MTG"
try:
    MT = trees.Tree("sample_mtg_forest.txt")
except trees.Tree, BT:
    MT = BT
    print MT.Display()
else:
    print "Failed to raise exception for reading a forest MTG"
# copy constructor of Tree
print "A copy of tree T:"
T2 = trees.Tree(MT)
print T2
print "Attribute names of copied tree", T2.Attributes()
