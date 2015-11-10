# -*- coding: utf-8 -*-
"""some tests for the class etrees.Trees: constructor and basic methods"""
__version__ = ""

import sys
import os
import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.trees.etrees as etrees
import openalea.stat_tool as stat_tool
from nose import with_setup


def init():
    return setup_func()

def setup_func():
    global  mtg_name

    mtg_name = "data/sample_mtg_forest.mtg"

    return mtg_name


@with_setup(setup_func)

def test_mtg_build():
    """constructor from a MTG"""
    msg1 = "Bad number of trees: "
    msg2 = "Bad number of variables: "
    mtg_t = etrees.Trees(mtg_name)
    msg1 += str(mtg_t.NbTrees())
    msg1 += " - should be 3"
    msg2 += str(mtg_t.NbVariables())
    msg2 += " - should be 4"
    assert mtg_t
    assert (mtg_t.NbTrees() == 3), msg1
    assert (mtg_t.NbVariables() == 4), msg2
    return mtg_t

def test_copy_vid_dic_list():
    """Copy the dictionaries of MTG and Tree vertex ids"""
    msg = "Copy error in the dictionaries of MTG and Tree vertex ids"
    T = test_mtg_build()
    TreeVidDictList = []
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictList.append(T.TreeVertexId(t))
        TreeVidDictListCp.append(T.TreeVertexId(t))
    T._SetMTGVidDictionary(TreeVidDictListCp)
    assert (T.TreeVertexId(0) == TreeVidDictList[0]), msg
    return T, TreeVidDictList

def test_copy_vid_bad_dic_list_no_check():
    """Copy the dictionaries of MTG and Tree vertex ids"""
    msg = "Copy error in the dictionaries of MTG and Tree vertex ids"
    T, TreeVidDictList = test_copy_vid_dic_list()
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictListCp.append(T.TreeVertexId(t))
        for k in TreeVidDictListCp[t].keys():
            TreeVidDictListCp[t][k] += 1
    T._SetMTGVidDictionary(TreeVidDictListCp)
    assert (T.TreeVertexId(0) == TreeVidDictListCp[0]), msg

def test_copy_vid_bad_dic_list_check():
    """Copy the dictionaries of MTG and Tree vertex ids"""
    msg = "Failed to raise IndexError"
    T, TreeVidDictList = test_copy_vid_dic_list()
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictListCp.append(T.TreeVertexId(t))
        for k in TreeVidDictListCp[t].keys():
            TreeVidDictListCp[t][k] += 1
    try:
        T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
    except IndexError, e:
        print e
    else:
        assert False, msg

def test_copy_vid_bad_value_dic_list_check():
    """Copy the dictionaries of MTG and Tree vertex ids.
    Check ValueError exception"""
    msg = "Failed to raise ValueError"
    T, TreeVidDictList = test_copy_vid_dic_list()
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictListCp.append(T.TreeVertexId(t))
    TreeVidDictListCp[0][3] = 0
    try:
        T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
    except ValueError, e:
        print e
    else:
        assert False, msg

def test_copy_vid_bad_type_dic_list_check():
    """Copy the dictionaries of MTG and Tree vertex ids.
    Check TypeError exception"""
    msg = "Failed to raise TypeError"
    T, TreeVidDictList = test_copy_vid_dic_list()
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictListCp.append(T.TreeVertexId(t))
        for k in TreeVidDictListCp[t].keys():
            TreeVidDictListCp[t][k] = str(TreeVidDictListCp[t][k])
    try:
        T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
    except TypeError, e:
        print e
    else:
        assert False, msg
    TreeVidDictListCp = []
    for t in range(T.NbTrees()):
        TreeVidDictListCp.append({})
        for k in TreeVidDictList[t].keys():
            TreeVidDictListCp[t][str(k)] = TreeVidDictList[t][k]
    try :
        T._SetMTGVidDictionary(TreeVidDictListCp, ValidityCheck=True)
    except Exception, msg:
        print msg
    else:
        assert False, msg

if __name__ == "__main__":
    mtg_name = init()
    test_mtg_build()
    test_copy_vid_bad_dic_list_no_check()
    test_copy_vid_bad_dic_list_check()
    test_copy_vid_bad_value_dic_list_check()
    test_copy_vid_bad_type_dic_list_check()