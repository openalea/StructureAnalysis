# -*- coding: utf-8 -*-
# create a MarkovOutTreeData with types
# {Short / Long} x {Vegetative / Floral}
import sys, os, datetime
import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt as hmt
import openalea.stat_tool as stat_tool
from stat_tool.plot import gnuplot, mplotlib
# stat_tool.plot.PLOTTER = mplotlib()
stat_tool.plot.PLOTTER = gnuplot()

T = trees.Trees("./Data/mtg_fuji_42_47.txt")
# 2 apple trees at annual shoot scale
# variable 0 : number of metamers
# variable 1 : state (0: veg. long 1: flo. long 2: veg. short 3: fl. short
# variable 2 : number of growth cycles
# variable 3 : is descendance unknown / censored? (case in 1999)

# get the dictionary of vertices with censored descencance
censored = []
for tid in range(T.NbTrees()):
    censored += [{}]
    TA = T.Tree(tid)
    # identifiers
    vids = TA.TreeVertex().values()
    for v in vids:
        # observed variables
        val = TA.Get(v)
        censored[tid][v] = (val[3] > 0)

def LengthFunc(x, TS=T):
    """Length function in PlantFrame"""
    return(TS.__mtg.node(x).NbInternode)
    
def ColorFunc(x, TS=T):
    """Color function in PlantFrame"""
    return(int(TS.__mtg.node(x).Type))

DiameterFunc = lambda x: 1

T.Plot("SojournSize", Title="Type", variable=1)

# T.MPlot("Counting", Title="Type", variable=1)

T.MPlot("Data", LengthFunc, ColorFunc, DiameterFunc, ColorFunc)

T.Plot("Data", LengthFunc, ColorFunc, DiameterFunc, ColorFunc)

# Keep variables 0 and 2 only
T2 = T.SelectVariable([0, 2], mode="Keep")
# Variable 1 should begin at value 0
T2 = T2.Shift(1,-1)

H = T2.Estimate("HIDDEN_MARKOV_TREE", 4, "IRREDUCIBLE", 0.9, 300)

