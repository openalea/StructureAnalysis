# -*- coding: utf-8 -*-
"""some tests for the class hmt.HiddenMarkovIndOutTree: constructor and basic methods"""
__version__ = ""

import sys
import os
# import openalea.tree_statistic.trees as trees
import openalea.tree_statistic.hmt as hmt
import openalea.stat_tool as stat_tool
from stat_tool.plot import gnuplot, mplotlib, DISABLE_PLOT

from nose import with_setup

def init():
    return setup_func()

def setup_func():
    stat_tool.plot.PLOTTER = mplotlib()
    stat_tool.plot.DISABLE_PLOT = True
    global hmt_name1, hmt_name2
    hmt_name1 = "data/hmot.hmt"
    hmt_name2 = "data/hmot_ascii.hmt"
    return hmt_name1, hmt_name2

@with_setup(setup_func)
def test_hmt_read():
    """Read some HMT from a file"""
    msg = "Could not read HMT file"
    try:
        H = hmt.HiddenMarkovIndOutTree(hmt_name1)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert H
        return H

def test_hmt_ascii_read():
    """Read some HMT from a file"""
    msg = "Could not read HMT file"
    try:
        H = hmt.HiddenMarkovIndOutTree(hmt_name2)
    except hmt.StatTreeError:
        assert False, msg
    else:
        assert H
        return H

def test_hmt_read_failure():
    """Use bad file name reading some HMT from a file"""
    msg = "Failed to raise exception for bad file name"
    try:
        H = hmt.HiddenMarkovIndOutTree("too_few_states_param.hmt")
    except hmt.StatTreeError, e:
        print e
    else:
        assert False, msg

def test_hmt_too_few_states_param_failure():
    """Read HMT from a file with too few states (parametric case)"""
    msg = "Failed to raise exception for too few state in parametric HMT"
    try:
        H = hmt.HiddenMarkovIndOutTree("")
    except hmt.StatTreeError, e:
        print e
    else:
        assert False, msg

def test_hmt_save():
    """Save some HMT into a file"""
    msg = "Could not save HMT into a file"
    H = test_hmt_read()
    file_name = "hmot_write.hmt"
    H.Save(file_name, "ASCII", True)
    try:
        f = file(file_name)
    except IOError:
        assert False, msg
    else:
        f.close()
        import os
        os.remove(file_name)
        return H

def test_hmt_spreadsheet_save():
    """Save some HMT into a file under spreadsheet formating"""
    msg = "Could not save HMT into a file"
    H = test_hmt_read()
    file_name = "hmot_spreadsheet_write.hmt"
    H.Save(file_name, "Spreadsheet", True)
    try:
        f = file(file_name)
    except IOError:
        assert False, msg
    else:
        f.close()
        import os
        os.remove(file_name)
        return H

def test_compare_display():
    """Compare display for two copies of an HMT"""
    msg = "Differences found in displaying equal HMTs: "
    equal = True
    differences = []
    H = test_hmt_read()
    H2 = hmt.HiddenMarkovIndOutTree(H)
    if str(H)!=str(H2):
        equal = False
    assert equal, msg + str(H) + "\n and" + str(H2)

def test_hmt_SAVEoREAD():
    """Compare Read and Save HMT"""
    msg = "load o (save o load) != load"
    equal = True
    H = test_hmt_read()
    file_name = "hmot_write.hmt"
    H.Save(file_name, "ASCII", True)
    try:
        f = file(file_name)
    except IOError:
        assert False, msg
    else:
        f.close()
    H2 = hmt.HiddenMarkovIndOutTree(file_name)
    if str(H)!=str(H2):
        equal = False
    import os
    os.remove(file_name)
    assert equal, msg

def test_extract_observation_histogram():
    H = test_hmt_read()
    HV = H.Extract("Observation", 1, 0)
    assert HV

def test_plot_first_occurrencer_histogram():
    H = test_hmt_read()
    if not(stat_tool.plot.DISABLE_PLOT):
        H.Plot(ViewPoint="FirstOccurrenceRoot", variable=1)

def test_plot_observation():
    H = test_hmt_read()
    if not(stat_tool.plot.DISABLE_PLOT):
        H.Plot()

if __name__ == "__main__":
    hmt_name1, hmt_name2 = init()
    stat_tool.plot.DISABLE_PLOT = False
    test_hmt_read()
    test_hmt_ascii_read
    test_hmt_read_failure()
    test_hmt_too_few_states_param_failure()
    test_hmt_save()
    test_hmt_spreadsheet_save()
    test_compare_display()
    test_hmt_SAVEoREAD()
    test_extract_observation_histogram()
    test_plot_first_occurrencer_histogram()
    test_plot_observation()
    # test_plot_data()
