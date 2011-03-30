# -*- coding: utf-8 -*-
"""some tests for the class trees.Trees: operations on trees"""
__version__ = ""

import sys
import os
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tool
from nose import with_setup


def init():
    return setup_func()

def setup_func():
    global  T,tv,nb_trees, tree_list, mtg_name
    # build some random initial tree
    stat_tool.plot.DISABLE_PLOT = True
    inf_bound = 0
    sup_bound = 3
    distrib = stat_tool.Uniform(inf_bound, sup_bound)
    max_depth = 3
    max_size = 10
    nbtrees = 40
    # define a set of trees
    tree_list = []
    tv = [1., 0, 1, 2.] # trees.TreeValue([1., 0])
    R = trees.TreeStructure(distrib, max_size, max_depth)
    tmp_tree = trees.Tree(tv, R)
    n = 1
    tree_list.append(trees.Tree(tmp_tree))
    while n < nbtrees:
        n = n+1
        R.Simulate(distrib, max_size, max_depth)
        tmp_tree = trees.Tree(tv, R)
        tree_list.append(trees.Tree(tmp_tree))
    distrib_list = []
    for i in range(tmp_tree.NbInt()):
        distrib_list.append(distrib)
    for n in range(len(tree_list)):
        tree_list[n].Simulate(distrib_list)
    T = trees.Trees(tree_list)
    nb_trees = nbtrees
    mtg_name = "data/sample_mtg_forest.mtg"
    
    return T, tv, nb_trees, tree_list, mtg_name


@with_setup(setup_func)

def test_mtg_build():
    """constructor from a MTG"""
    mtg_t = trees.Trees(mtg_name)
    assert mtg_t
    return mtg_t

def test_cluster():
    """cluster values of variable 0"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.Cluster("Step", 0, 10)
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)
    clust = T

def test_cluster_step_failure():
    """use invalid step in clustering values"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.Cluster("Step", 0, -10)
    except trees.StatTreeError, msg:
        print msg
    else:
        msg = "Failed to raise exception for step in clustering"
        assert False, msg

def test_difference():
    """1st order differenciation for variable 0"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.Difference(0)
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)

def test_merge():
    """merge a tree with itself"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.Merge([mtg_t])
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)
    msg1 = "Bad number of merged trees"
    assert T.NbTrees() == mtg_t.NbTrees()*2, msg1

def test_merge_variables():
    """merge variables of a tree with themselves"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.MergeVariable([mtg_t])
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)
    msg1 = "Bad number of Variables"
    assert (T.NbVariables() == mtg_t.NbVariables()*2), msg1
    msg2 = "Error copying correspondence between Tree and MTG vertex ids: "
    msg2 += str(T.MTGVertexId(0)) + " - should be "
    msg2 += str(mtg_t.MTGVertexId(0))
    assert (mtg_t.MTGVertexId(0) == T.MTGVertexId(0)), msg2

def test_select_variable():
    """select variable 0"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.SelectVariable(0)
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)
    msg = "Bad number of trees after variable selection"
    assert T.NbTrees() == mtg_t.NbTrees(), msg

def test_merge_failure():
    """use invalid tree in merge"""
    mtg_t = test_mtg_build()
    try:
        TS = mtg_t.SelectVariable(0)
        T = mtg_t.Merge([TS])
    except trees.StatTreeError, msg:
        print msg
    else:
        msg = "Failed to raise exception for bad number "+\
            "of variables in merge"
        assert False, msg

def test_build_vectors():
    vec = T.BuildVectors()
    assert vec
    import openalea
    assert type(vec) == openalea.stat_tool._stat_tool._Vectors

if __name__ == "__main__":
    T, tv, nb_trees, tree_list, mtg_name = init()
    test_mtg_build()
    test_cluster()
    test_cluster_step_failure()
    test_difference()
    test_merge()
    test_select_variable()
    test_merge_failure()
    test_build_vectors()
