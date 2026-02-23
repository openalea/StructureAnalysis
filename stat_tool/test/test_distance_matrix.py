"""distance matrix tests


:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""

__version__ = "$Id$"

from .tools import robust_path as get_shared_data
from .tools import interface

from openalea.stat_tool.cluster import Clustering, ToDistanceMatrix
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.vectors import VectorDistance, Vectors

import pytest

@pytest.fixture
def data():
    # Inspired from stat_toot_test.aml
    vec10 = Vectors(get_shared_data("chene_sessile.vec"))
    vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
    matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
    c1 = Clustering(
        matrix10, "Partition", 3, Prototypes=[1, 3, 12], Algorithm="Divisive"
        )
    return ToDistanceMatrix(c1)

@pytest.fixture
def myi(data):
    return interface(data, "data/distribution1.dist", ToDistanceMatrix)


def test_print(myi):
    myi.print_data()


def test_plot(myi):
    myi.plot()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()


def test_symmetrize(data):
    assert data.symmetrize()

def test_unnormalize(data):
    assert data.symmetrize().unnormalize()

def test_select_individual(data):
    keep_false = data.select_individual([1], keep=False)
    keep_true = data.select_individual([1], keep=True)

    assert keep_false.nb_row == 2
    assert keep_false.nb_column == 2

    assert keep_true.nb_row == 1
    assert keep_true.nb_column == 1

def test_get_distance(data):
    assert data.get_distance(0, 0)

def test_get_length(data):
    assert data.get_length(0, 0)

def test_get_substitution(data):
    # no substituion computed for vectors
    assert data.get_substitution_distance(0, 0) == -1

def test_get_insertion(data):
    # no insertion computed for vectors
    assert data.get_insertion_distance(0, 0) == -1

def test_get_transposition(data):
    # no transposition computed for vectors
    assert data.get_transposition_distance(0, 0) == -1

def test_get_deletion(data):
    # no deletion computed for vectors
    assert data.get_deletion_distance(0, 0) == -1

def test_get_length_outside_range(data):
    try:
        data.get_length(-1, 0)
        assert False
    except TypeError:
        assert True

