# -*- coding: utf-8 -*-
"""some tests for the class hmt.HiddenMarkovIndOutTree: basic operations"""
__version__ = ""

import sys
import os
# import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt as hmt
import openalea.stat_tool as stat_tool
from nose import with_setup

def init():
    return setup_func()

def setup_func():
    global hmt_name, H
    hmt_name = "data/hmot_cycle.hmt"
    H = hmt.HiddenMarkovIndOutTree(hmt_name)
    return hmt_name, H

@with_setup(setup_func)
def test_nb_states():
    """Test the number of states"""
    msg = "Bad number of states"
    assert H.NbStates() == 5, msg

def test_hmt_data_failure():
    """Extract data part from plain HMT"""
    msg = "Failed to raise exception for missing data part"
    try:
        H.ExtractData()
    except hmt.StatTreeError, e:
        print e
    else:
        assert False, msg

def build_binary_tree():
    """Build a binary tree"""
    distrib = stat_tool.Uniform(2, 2)
    max_depth = 5
    max_size = 200
    nbtrees = 1
    # define a set of trees
    tv = [0]
    import openalea.tree_statistic.trees as trees
    R = trees.TreeStructure(distrib, max_size, max_depth)
    tmp_tree = trees.Tree(tv, R)
    T = trees.Trees([tmp_tree])
    return T

def test_hmt_simulate():
    """Simulate a forest from HMT"""
    msg1 = "Bad value for initial state"
    msg2 = "Bad value for final observation"
    msg3 = "Bad value observation at state 1"
    T = build_binary_tree()
    THMT = H.Simulate(T)
    assert(THMT.Tree(0).Get(0)[0]==0), msg1
    assert(THMT.Tree(0).Get(31)[1]==0), msg2
    assert(THMT.Tree(0).Get(7)[1] in [6, 7]), msg3
    return T, THMT
    
def test_extract_model():
    """Extract model part of HiddenMarkovTreeData"""
    T, THMT = test_hmt_simulate()
    # THMT.Display(Detail=2)
    H2 = THMT.ExtractMarkov()
    assert H2, msg
    THMT2 = H2.ExtractData()
    msg = "Bad value for extracted HiddenMarkovTree"
    file_name = "hmot_extract_write.hmt"
    H2.Save(file_name, "ASCII", True)
    H3 = hmt.HiddenMarkovIndOutTree(hmt_name)
    import os
    os.remove(file_name)
    s1 = H3._chmt().Display(True)
    s2 = H._chmt().Display(True)
    assert(s1 == s2), msg
    assert(str(THMT)==str(THMT2)), msg

def test_hmt_permutation():
    """Permute the states of an HMT"""
    P = [1, 2, 3, 4, 0]
    # state i is renamed i+1,
    # state 4 is renamed 0 
    H2 = hmt.HiddenMarkovIndOutTree(H)
    H2.StatePermutation(P)
    T = build_binary_tree()
    THMT = H2.Simulate(T)
    msg1 = "Bad value for initial state"
    msg2 = "Bad value for final observation"
    msg3 = "Bad value observation at state 1"
    assert(THMT.Tree(0).Get(0)[0]==1), msg1
    assert(THMT.Tree(0).Get(31)[1]==0), msg2
    assert(THMT.Tree(0).Get(7)[1] in [6, 7]), msg3
    return T, THMT

if __name__ == "__main__":
    hmt_name, H = init()
    test_nb_states()
    test_hmt_data_failure()
    test_hmt_simulate()
    test_extract_model()
    test_hmt_permutation()
    
