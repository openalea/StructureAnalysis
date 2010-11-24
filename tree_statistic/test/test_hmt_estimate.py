# -*- coding: utf-8 -*-
"""some tests for the class hmt.HiddenMarkovTree: estimation"""
__version__ = ""

import sys
import os
import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt as hmt
import openalea.stat_tool as stat_tool
from nose import with_setup

def init():
    return setup_func()

def setup_func():
    global hmt_name, hmt_name, T, HInit
    mtg_name = "data/sim3.mtg"
    hmt_name = "data/hmot3_init.hmt"
    HInit = hmt.HiddenMarkovTree(hmt_name)
    T = trees.Trees(mtg_name)
    return hmt_name, mtg_name, T, HInit

@with_setup(setup_func)
def test_estimate_init_hmt():
    """Estimate with an initial HMT read from a file"""
    msg = "Could not estimate HiddenMarkovTree"
    try:        
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_state():
    """Estimate with an initial HMT specified by initial transitions"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", 3, "IRREDUCIBLE", 0.9, 50)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_hmt_options():
    """Estimate with an initial HMT read from a file, using options"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50, "ForwardBackward", True)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_hmt_value():
    """Estimate with an initial HMT read from a file, and compare result
    with a save HMT model"""
    msg = "Bad estimation for HiddenMarkovTree"
    hmtref_name = "data/hmot3estim.hmt"
    HREF = hmt.HiddenMarkovTree(hmtref_name)
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50, "ForwardBackward", True)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert (EH.Display() == HREF.Display())
        return EH

def test_estimate_init_state_options():
    """Estimate with an initial HMT pecified by initial transitions, using options"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", 3, "IRREDUCIBLE", 0.9, 50, "ForwardBackward", False)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_hmt_saem():
    """Estimate with an initial HMT read from a file, using SAEM"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50, "ForwardBackward", True,
                        Algorithm="ForwardBackwardSampling", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_state_gibbs():
    """Estimate with an initial HMT pecified by initial transitions, using gibbs sampling"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", 3, "IRREDUCIBLE", 0.9, 50, "ForwardBackward", False,
                        Algorithm="GibbsSampling", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_hmt_gibbs():
    """Estimate with an initial HMT read from a file, using SAEM"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50, "ForwardBackward", True,
                        Algorithm="GibbsSampling", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH


def test_estimate_init_state_cem():
    """Estimate with an initial HMT pecified by initial transitions, using Viterbi restoration"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", 3, "IRREDUCIBLE", 0.9, 50, "ForwardBackward", False,
                        Algorithm="Viterbi", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_hmt_cem():
    """Estimate with an initial HMT read from a file, using SAEM"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", HInit, 50, "ForwardBackward", True,
                        Algorithm="Viterbi", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_estimate_init_state_saem():
    """Estimate with an initial HMT pecified by initial transitions, using SAEM"""
    msg = "Could not estimate HiddenMarkovTree"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", 3, "IRREDUCIBLE", 0.9, 50, "ForwardBackward", False,
                        Algorithm="ForwardBackwardSampling", Saem=0.5)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert EH
        return EH

def test_bad_hmt_init_failure():
    """Estimate HMT with initial HMT of incorrect type"""
    msg = "Failed to raise exception for incorrect initial HMT"
    try:
        EH = T.Estimate("HIDDEN_MARKOV_TREE", hmt, 50)
    except TypeError, e:
        print e
    else:
        assert False, msg


if __name__ == "__main__":
    hmt_name, mtg_name, T, HInit = init()
    test_estimate_init_hmt()
    test_estimate_init_state()
    test_estimate_init_hmt_options()
    test_estimate_init_state_options()
    test_estimate_init_hmt_saem()
    test_estimate_init_state_saem()
    test_estimate_init_hmt_gibbs()
    test_estimate_init_state_gibbs()
    test_estimate_init_hmt_cem()
    test_estimate_init_state_cem()
    test_bad_hmt_init_failure()
