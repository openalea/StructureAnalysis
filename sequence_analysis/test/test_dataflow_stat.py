from openalea.core.alea import *
# !!important!! import dataflowview, which defines the fields of each nodes
#from openalea.grapheditor import dataflowview

pm = PackageManager()
pm.init(verbose=True)

# These tests use gnuplot interface, which requires human interaction
# Consequently, they cannot be used within builbot (which hangs forever)
# We added a flags inside aml/src/aml/wralea/py_stat.py to prevent gnuplot
# to be launched if these tests are run with nosetests. The remaining of the
# nodes are run.
# In order to have the gnuplot interface, run this script with python instead of nosetests

def test_demo_corsican():
    """ Test changepoint demo corsican  """
    res = run(('demo.changepoint_stat_tool','Corsican pine change point'),{},pm=pm)
    assert res == []

def test_demo_dycorinia():
    """ Test dataflow demo dycorinia """
    res = run(('demo.changepoint_stat_tool','Dycorinia change point'),{},pm=pm)
    assert res == []

def test_oak_demo():
    """ Test dataflow demo oak"""
    res = run(('demo.changepoint_stat_tool', 'oak_demo'),{},pm=pm)
    assert res == []

def test_beech1_demo():
    """ Test dataflow demo beech"""
    res = run(('demo.changepoint_stat_tool', 'beech1'),{},pm=pm)
    assert res == []

def test_stat_tool_demos_and_tutorial_convolution():
    """ Test dataflow demo compound tutorial"""
    res = run(('demo.stat_tool demos and tutorials', 'convolution_tutorial'),{},pm=pm)
    assert res == []


def test_stat_tool_demos_and_tutorial_compound():
    res = run(('demo.stat_tool demos and tutorials', 'compound_tutorial'),{},pm=pm)
    assert res == []



if __name__ == "__main__":
    test_demo_corsican()
    test_demo_dycorinia()
    test_stat_tool_tutorial_compound()
    test_stat_tool_tutorial_convolution()

