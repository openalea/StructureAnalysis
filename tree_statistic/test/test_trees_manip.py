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
    global nb_trees, tree_list, tv, T, mtg_name
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
    # return nb_trees, tree_list, tv, T, mtg_name
    return mtg_name
    
def check_mtg_mapping(mtg_t):
    """Check consistency between tree vertices
    and mapping betweeen mtg vertices and tree vertices"""
    for t in range(mtg_t.NbTrees()):
        # map: TreeId -> MTGVertexId
        Tr = mtg_t.Tree(t)
        tmap = Tr.MTGVertex()
        tmapinv = Tr.TreeVertex()
        for m in [tmap, tmapinv]:
            if len(m) != Tr.Size():
                msg = "Length of mapping (" + str(len(m))
                msg += ") and number of tree vertices (" + str(Tr.Size())
                msg += ") do not match for tree: " + str(t)
                assert len(m) == Tr.Size(), msg
        for i in tmap.keys():
            if tmapinv[tmap[i]] != i:
                msg = "Inconsistency in mapping and inverse mapping"
                msg += " in tree: " + str(t)
                assert tmapinv[tmap[i]] != i, msg

@with_setup(setup_func)

def test_mtg_build():
    """constructor from a MTG"""
    mtg_t = trees.Trees(mtg_name)
    assert mtg_t
    check_mtg_mapping(mtg_t)
    return mtg_t

def test_tree_size():
    """test size of trees"""
    msg = "Bad total number of vertices: "
    mtg_t = test_mtg_build()
    ts = mtg_t.Size()
    msg += str(ts)
    ts2 = 0
    for i in range(mtg_t.NbTrees()):
        ts2 += mtg_t.Size(i)

    msg += " - should be: " + str(ts2)
    assert (ts == ts2), msg
    assert (ts == 43), msg

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
    msg = "Bad number of variables after variable selection"
    assert T.NbVariables() == 1, msg

def test_select_individual():
    """Select individuals 0, 2"""
    mtg_t = test_mtg_build()
    try:
        T = mtg_t.SelectIndividual([0, 2])
    except trees.StatTreeError, msg:
        assert False, str(msg) + "\n" + str(mtg_t)
    msg = "Bad number of trees after variable selection"
    assert T.NbTrees() == 2, msg
    msg = "Bad number of variables after variable selection"
    assert T.NbVariables() == mtg_t.NbVariables(), msg
    check_mtg_mapping(T)
    return T

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
    """Extract vectors from a tree"""
    msg = "Bad number of vectors: "
    vec = T.BuildVectors()
    assert vec
    import openalea
    assert type(vec) == openalea.stat_tool._stat_tool._Vectors
    msg += str(vec.nb_vector)
    msg += " - should be :"
    msg += str(T.Size())
    assert (vec.nb_vector == T.Size()), msg

def test_build_sequences():
    """Extract non-redundant sequences from a tree.
       Cut trees at '+' edges"""
    mtg_t = test_mtg_build()
    from openalea.tree_statistic.trees import etrees
    for t in range(mtg_t.NbTrees()):
        msg = "Bad number of items in tree " + str(t)
        msg += ": "
        Tr = etrees.Tree(mtg_t.Tree(t))
        T = etrees.Trees([Tr])
        seq = T.BuildPySequences(False, False)
        assert seq
        nb_items = seq.build_vectors(True).nb_vector
        msg += str(nb_items)
        msg += " - should be :"
        msg += str(T.Size())
        assert (nb_items == T.Size()), msg
        e = [v for v in Tr.Postorder() if not(Tr.IsRoot(v)) \
            and Tr.EdgeType(Tr.Parent(v), v) == "+"]
        nb_seq = len(e) + 1
        msg = "Bad number of sequences in tree " + str(t)
        msg += ": " + str(len(seq))
        msg += " - should be: " + str(nb_seq)
        assert len(seq) == nb_seq, msg

def test_build_auto_axis_sequences():
    """Extract non-redundant sequences from a tree.
       Cut trees at '+' edges. Unique children are
       considered as '<' in every case."""
    mtg_t = test_mtg_build()
    from openalea.tree_statistic.trees import etrees
    for t in range(mtg_t.NbTrees()):
        msg = "Bad number of items in tree " + str(t)
        msg += ": "
        Tr = etrees.Tree(mtg_t.Tree(t))
        T = etrees.Trees([Tr])
        seq = T.BuildPySequences(False, True)
        assert seq
        nb_items = seq.build_vectors(True).nb_vector
        msg += str(nb_items)
        msg += " - should be :"
        msg += str(T.Size())
        assert (nb_items == T.Size()), msg
        e = [v for v in Tr.Postorder() if not(Tr.IsRoot(v)) and \
            (Tr.EdgeType(Tr.Parent(v), v) == "+") and \
            (Tr.NbChildren(Tr.Parent(v)) > 1)]
        nb_seq = len(e) + 1
        msg = "Bad number of sequences in tree " + str(t)
        msg += ": " + str(len(seq))
        msg += " - should be: " + str(nb_seq)
        assert len(seq) == nb_seq, msg

def test_build_redundant_sequences():
    """Extract redundant sequences from a tree.
       Follow every path from a leaf to root vertex."""
    mtg_t = test_mtg_build()
    from openalea.tree_statistic.trees import etrees
    for t in range(mtg_t.NbTrees()):
        msg = "Bad number of items in tree " + str(t)
        msg += ": "
        Tr = etrees.Tree(mtg_t.Tree(t))
        T = etrees.Trees([Tr])
        seq = T.BuildPySequences(True)
        assert seq
        e = [v for v in Tr.Postorder() if 
             Tr.NbChildren(v) == 0]
        nb_seq = len(e)
        msg = "Bad number of sequences in tree " + str(t)
        msg += ": " + str(len(seq))
        msg += " - should be: " + str(nb_seq)
        assert len(seq) == nb_seq, msg
        
def test_extract_mtg():
    """Extract an MTG object from a trees"""
    msg1 = "Extraction of MTG failed"
    msg2 = "Bad number of trees: "
    g = trees.Mtg(T)
    assert g, msg1
    assert len(g.vertices(scale=1)) == T.NbTrees(), msg2
    
if __name__ == "__main__":
    # nb_trees, tree_list, tv, T, mtg_name = init()
    mtg_name = init()
    test_mtg_build()
    test_tree_size()
    test_cluster()
    test_cluster_step_failure()
    test_difference()
    test_merge()
    test_select_variable()
    test_select_individual()
    test_merge_failure()
    test_build_vectors()
    test_build_sequences()
    test_build_auto_axis_sequences()
    test_build_redundant_sequences()
    test_extract_mtg()
