# -*- coding: utf-8 -*-
"""some tests for the class trees.Trees: constructor and basic methods"""
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
def test_types():
    """Test the variables types"""
    msg = "Type lengths do not match"
    print(T.Types())
    assert len(T.Types()) == len(tv), msg

def test_display():
    print(T.Display(Detail=1))
    print(T.Display(Detail=2))
    T.Display(ViewPoint="DATA", Detail=1)

def test_tree():
    """Extraction of trees"""
    msg = "Error extracting tree 0"
    equal = True
    differences = []
    if str(T.Tree(0))!=str(tree_list[0]):
        equal = False
    assert equal, msg

def test_compare_display():
    """Compare the trees contained in T
        and those used for builing it"""
    msg = "Differences found in displaying equal trees: "
    equal = True
    differences = []
    for n in range(len(tree_list)):
        if str(T.Tree(n))!=str(tree_list[n]):
            equal = False
            differences.append(n)
    assert equal, msg + str(differences)

def test_extract_histogram():
    HV = T.ExtractHistogram("Value", variable=1)
    assert HV

def test_plot_histogram():
    T.Plot(ViewPoint="FirstOccurrenceRoot", variable=1)

def test_plot_data():
    if not(stat_tool.plot.DISABLE_PLOT):
        T.Plot()

def test_nb_trees():
    """test for the number of trees"""
    msg = "Bad number of trees"
    assert nb_trees == T.NbTrees(), msg

def test_copy():
    """test the copy of trees"""
    T2 = trees.Trees(T)
    msg = "Copy error"
    assert T2._ctrees_display() == T._ctrees_display(), msg

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
    except trees.FormatError, msg:
        assert False, str(msg) + "\n" + str(T.Tree(0))
    clust = T

def test_build_vectors():
    vec = T.BuildVectors()
    assert vec
    import openalea
    assert type(vec) == openalea.stat_tool._stat_tool._Vectors

def test_build_mtg_filter():
    """Read a Tree from a MTG with filter and custom attributes"""
    import openalea.aml as amlPy
    f = lambda x: amlPy.Feature(x, "Length")*3*(amlPy.Feature(x, "Diam")/2)**2
    filter = lambda x: x < 6
    attributes = ["anything"]
    T = trees.Trees(mtg_name, filter, attributes, [f], scale=2)
    msg = "Error building tree from MTG using filters"
    assert T, msg
    assert T.NbTrees() == 1

def test_mtg_build_failure():
    """use invalid MTG file name"""
    try:
        T = trees.Trees("no_such_file.t\/t")
    except IOError, e:
        print e
    else:
        msg = "Failed to raise exception for invalid MTG file name"
        assert False, msg

def test_build_attribute_failure():
    """use inconsistent attribute number in building trees"""
    f = lambda x: 1
    filter = lambda x: True
    try:
        T = trees.Trees(mtg_name, filter, [], [f], scale=2)
    except ValueError, v:
        print v
    else:
        msg = "Failed to raise exception for inconsistent attribute number"
        assert False, msg

def test_attribute_name_failure():
    """use bad attribute name in building trees"""
    """does not seem to be run after all"""
    f = lambda x: 1
    filter = lambda x: True
    try:
        T = trees.Trees(mtg_name, filter, [f], [f], scale=2)
    except TypeError, t:
        print t
    else:
        msg = "Failed to raise exception for bad attribute name"
        assert False, msg

def tst_attribute_function_failure():
    """use bad attribute function in building trees"""
    filter = lambda x: True
    attributes = ["anything"]
    try:
        T = trees.Trees(mtg_name, filter, attributes, attributes,
                        scale=2)
    except TypeError, t:
        print t
    else:
        msg = "Failed to raise exception for bad attribute function"
        assert False, msg

def test_attribute_type_failure():
    """use bad attribute types in building trees"sample_mtg_forest.txt"""
    filter = lambda x: x < 6
    attributes = ["anything"]
    import openalea.aml as amlPy
    try:
        T = trees.Trees(mtg_name, filter, attributes,
                        [lambda x: "z"], scale=2)
    except TypeError, t:
        print t
    else:
        print "Failed to raise exception for bad attribute type"
    return 0

def test_attribute_name_failure():
    """use bad attribute names in building trees"""
    filter = lambda x: x < 6
    f = lambda x: 1
    try:
        T = trees.Trees(mtg_name, filter, [f], [f], scale=2)
    except TypeError, t:
        print t
    else:
        msg = "Failed to raise exception for bad attribute name"
        assert False, msg
    return 0

def test_filter_failure():
    """use bad filterin function in building trees"""
    filter = lambda x: x < 6
    f = lambda x: 1
    attributes = ["anything"]
    try:
        T = trees.Trees(mtg_name, attributes, attributes, [f], scale=2)
    except TypeError, t:
        print t
    else:
        msg = "Failed to raise exception for bad attribute function"
        assert False, msg

def test_nonrecursive_filter_failure():
    """use nonrecursive filtering function in building trees"""
    filter = lambda x: x < 6
    f = lambda x: 1
    attributes = ["anything"]
    try:
        T = trees.Trees(mtg_name, lambda x: x != 2, attributes, [f],
                        scale=2)
    except IndexError, i:
        print i
    else:
        msg = "Failed to raise exception for filter not filtering descendants"
        assert False, msg

def test_filter_type_failure():
    """use incorrect type return for filter in building trees"""
    filter = lambda x: x < 6
    f = lambda x: "a"
    attributes = ["anything"]
    try:
        T = trees.Trees(mtg_name, filter, attributes, [f], scale=1)
    except TypeError, t:
        print t
    else:
        msg = "Failed to raise exception for bad filter type"
        assert False, msg

if __name__ == "__main__":
    T, tv, nb_trees, tree_list, mtg_name = init()
    test_types()
    test_display()
    test_tree()
    test_compare_display()
    test_extract_histogram()
    test_plot_histogram()
    test_plot_data()
    test_nb_trees()
    test_copy()
    test_mtg_build()
    test_cluster()
    test_mtg_build_failure()
    test_attribute_name_failure()
    test_attribute_function_failure()
    test_attribute_type_failure()
    test_attribute_name_failure()
    test_filter_failure()
    test_nonrecursive_filter_failure()
    test_filter_type_failure()
