"""Compound tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__version__ = "$Id$"

from .tools import interface, robust_path as get_shared_data

from openalea.stat_tool.compound import Compound
from openalea.stat_tool.data_transform import ExtractDistribution
from openalea.stat_tool.distribution import Binomial, NegativeBinomial, set_seed
from openalea.stat_tool.estimate import Estimate

import pytest

@pytest.fixture
def data():
    d1 = Binomial(2, 5, 0.5)
    d2 = NegativeBinomial(0, 2, 0.5)
    comp = Compound(d1, d2)
    return comp

@pytest.fixture
def myinterface(data):
    return interface(data, "data/compound1.cd", Compound)


def test_constructor_from_dists_and_threshold(data):
    compound1 = data
    compound2 = Compound(
        Binomial(2, 5, 0.5), NegativeBinomial(0, 2, 0.5), Threshold=0.99999
    )
    assert str(compound1) == str(compound2)

def test_print(myinterface):
    myinterface.print_data()

def test_display(myinterface):
    myinterface.display()
    myinterface.display_versus_ascii_write()
    myinterface.display_versus_str()


def test_save(myinterface):
    myinterface.save()

def test_plot_write(myinterface):
    myinterface.plot_write()

def test_file_ascii_write(myinterface):
    myinterface.file_ascii_write()

def test_spreadsheet_write(myinterface):
    myinterface.spreadsheet_write()

def test_simulate(myinterface):
    sim = myinterface.simulate()
    sim.plot()

def test_extract(data):
    """run and test the extract methods"""
    m = data
    assert m.extract_compound() == ExtractDistribution(m, "Compound")
    assert m.extract_sum() == Binomial(2, 5, 0.5)
    assert m.extract_sum() == ExtractDistribution(m, "Sum")
    l1 = str(m.extract_elementary())
    l2 = str(NegativeBinomial(0, 2, 0.5))
    assert str(l1.split("\n")[0:3]) == str(l2.split("\n")[0:3])
    assert m.extract_elementary() == ExtractDistribution(m, "Elementary")

def test_extract_data(myinterface):
    """todo : check if this test makes sense"""

    s = myinterface.simulate()
    # e = Estimate(s, "Compound",  Binomial(2, 5, 0.5), "Sum")
    d = s.extract_sum()
    assert d
    _eprime = Estimate(s, "COMPOUND", Binomial(0, 10, 0.5), "Sum")

def test_truncate(data):
    s = data
    _res = s.truncate(4)


