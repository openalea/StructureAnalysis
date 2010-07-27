import runpy


def _test_functional1():
    runpy.run_module('functional1')

def test_functional2():
    runpy.run_module('functional2')

def _test_functional3():
    runpy.run_module('functional3')
