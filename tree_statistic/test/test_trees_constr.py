# -*- coding: utf-8 -*-
"""some tests for the class trees.Trees: constructor and basic methods"""
__version__ = ""

import sys
import os
import openalea.tree_statistic.trees as trees
import openalea.stat_tool as stat_tool


class Test():
    """a simple unittest class
    """
    def __init__(self):
        """build some random initial tree"""
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
        self.T = trees.Trees(tree_list)
        self.tv = tv
        self.nb_trees = nbtrees
        self.tree_list = tree_list
        self.mtg_t = []
        self.clust = []
        self.mtg_name = "data/sample_mtg_forest.mtg"

    def test_types(self):
        """Test the variables types"""
        T = self.T
        msg = "Type lengths do not match"
        print(T.Types())
        assert len(T.Types()) == len(self.tv), msg
        return T.Types()

    def test_display(self):
        T = self.T
        print(T.Display(Detail=1))
        print(T.Display(Detail=2))
        T.Display(ViewPoint="DATA", Detail=1)

    def test_tree(self):
        """Extraction of trees"""
        tree_list = self.tree_list
        T = self.T
        msg = "Error extracting tree 0"
        equal = True
        differences = []
        if str(T.Tree(0))!=str(tree_list[0]):
            equal = False
        assert equal, msg

    def test_compare_display(self):
        """Compare the trees contained in self.T
            and those used for builing it"""
        tree_list = self.tree_list
        T = self.T
        msg = "Differences found in displaying equal trees: "
        equal = True
        differences = []
        for n in range(len(tree_list)):
            if str(T.Tree(n))!=str(tree_list[n]):
                equal = False
                differences.append(n)
        assert equal, msg + str(differences)

    def test_extract_histogram(self):
        T = self.T
        HV = T.ExtractHistogram("Value", variable=1)
        assert HV

    def test_plot_histogram(self):
        self.T.Plot(ViewPoint="FirstOccurrenceRoot", variable=1)

    def test_plot_data(self):
        if not(stat_tool.plot.DISABLE_PLOT):
            self.T.Plot()

    def test_nb_trees(self):
        """test for the number of trees"""
        T = self.T
        msg = "Bad number of trees"
        assert self.nb_trees == T.NbTrees(), msg

    def test_copy(self):
        """test the copy of trees"""
        T = self.T
        T2 = trees.Trees(T)
        msg = "Copy error"
        assert T2._ctrees_display() == T._ctrees_display(), msg

    def test_mtg_build(self):
        """constructor from a MTG"""
        self.mtg_t = trees.Trees(self.mtg_name)
        assert self.mtg_t

    def test_cluster(self):
        """cluster values of variable 0"""
        T = self.mtg_t
        if type(T) == list:
            self.test_mtg_build()
        try:
            T = self.mtg_t.Cluster("Step", 0, 10)
        except trees.FormatError, msg:
            assert False, str(msg) + "\n" + str(T.Tree(0))
        self.clust = T

    def test_build_vectors(self):
        T = self.clust
        if type(T) == list:
            self.test_cluster()
        vec = self.T.BuildVectors()
        assert vec
        import openalea
        assert type(vec) == openalea.stat_tool._stat_tool._Vectors

    def test_build_mtg_filter(self):
        """Read a Tree from a MTG with filter and custom attributes"""
        import openalea.aml as amlPy
        if type(self.mtg_t) == list:
            self.test_mtg_build()
        f = lambda x: amlPy.Feature(x, "Length")*3*(amlPy.Feature(x, "Diam")/2)**2
        filter = lambda x: x < 6
        attributes = ["anything"]
        T = trees.Trees(self.mtg_name, filter, attributes, [f], scale=2)
        msg = "Error building tree from MTG using filters"
        assert T
        assert T.NbTrees() == 1

    def test_mtg_build_failure(self):
        """use invalid MTG file name"""
        try:
            T = trees.Trees("no_such_file.t\/t")
        except IOError, e:
            print e
        else:
            msg = "Failed to raise exception for invalid MTG file name"
            assert False, msg

    def test_build_attribute_failure(self):
        """use inconsistent attribute number in building trees"""
        f = lambda x: 1
        filter = lambda x: True
        try:
            T = trees.Trees(self.mtg_name, filter, [], [f], scale=2)
        except ValueError, v:
            print v
        else:
           msg = "Failed to raise exception for inconsistent attribute number"
           assert False, msg

    def test_attribute_name_failure(self):
        """use bad attribute name in building trees"""
        """does not seem to be run after all"""
        f = lambda x: 1
        filter = lambda x: True
        try:
            T = trees.Trees(self.mtg_name, filter, [f], [f], scale=2)
        except TypeError, t:
            print t
        else:
            msg = "Failed to raise exception for bad attribute name"
            assert False, msg

    def test_attribute_function_failure(self):
        """use bad attribute function in building trees"""
        #filter = lambda x: True
        #attributes = ["anything"]
        #try:
            #T = trees.Trees(self.mtg_name, filter, attributes, attributes,
                            #scale=2)
        #except TypeError, t:
            #print t
        #else:
            #msg = "Failed to raise exception for bad attribute function"
            #assert False, msg
        return 0

    def test_attribute_type_failure(self):
        """use bad attribute types in building trees"sample_mtg_forest.txt"""
        #filter = lambda x: x < 6
        #attributes = ["anything"]
        #if type(self.mtg_t) == list:
            #self.test_mtg_build()
        #import openalea.aml as amlPy
        #try:
            #T = trees.Trees(self.mtg_name, filter, attributes,
                            #[lambda x: amlPy.Descendants(x)], scale=2)
        #except TypeError, t:
            #print t
        #else:
            #print "Failed to raise exception for bad attribute type"
        return 0

    def test_attribute_name_failure(self):
        """use bad attribute names in building trees"""
        #filter = lambda x: x < 6
        #f = lambda x: 1
        #if type(self.mtg_t) == list:
            #self.test_mtg_build()
        #import openalea.aml as amlPy
        #try:
            #T = trees.Trees(self.mtg_name, filter, [f], [f], scale=2)
        #except TypeError, t:
            #print t
        #else:
            #msg = "Failed to raise exception for bad attribute name"
            #assert False, msg
        return 0

    def test_filter_failure(self):
        """use bad filterin function in building trees"""
        filter = lambda x: x < 6
        f = lambda x: 1
        if type(self.mtg_t) == list:
            self.test_mtg_build()
        import openalea.aml as amlPy
        attributes = ["anything"]
        try:
            T = trees.Trees(self.mtg_name, attributes, attributes, [f], scale=2)
        except TypeError, t:
            print t
        else:
            msg = "Failed to raise exception for bad attribute function"
            assert False, msg

    def test_nonrecursive_filter_failure(self):
        """use nonrecursive filtering function in building trees"""
        filter = lambda x: x < 6
        f = lambda x: 1
        if type(self.mtg_t) == list:
            self.test_mtg_build()
        import openalea.aml as amlPy
        attributes = ["anything"]
        try:
            T=trees.Trees(self.mtg_name, lambda x: x != 2, attributes, [f],
                          scale=2)
        except IndexError, i:
            print i
        else:
            msg = "Failed to raise exception for filter not filtering descendants"
            assert False, msg

    def test_filter_type_failure(self):
        """use incorrect type return for filter in building trees"""
        #filter = lambda x: x < 6
        #f = lambda x: 1
        #if type(self.mtg_t) == list:
            #self.test_mtg_build()
        #import openalea.aml as amlPy
        #attributes = ["anything"]
        #try:
            #T = trees.Trees(self.mtg_name, filter, attributes, [f], scale=1)
        #except TypeError, t:
            #print t
        #else:
            #msg = "Failed to raise exception for bad filter type"
            #assert False, msg
        return 0


if __name__ == "__main__":
    runTestClass(Test())

