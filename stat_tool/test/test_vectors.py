"""vectors tests"""

from .tools import interface, runTestClass
from .tools import robust_path as get_shared_data

from openalea.stat_tool.enums import variance_type
from openalea.stat_tool.vectors import (
    ComputeRankCorrelation,
    ContingencyTable,
    OutputFormat,
    VarianceAnalysis,
    Vectors,
)

import pytest

@pytest.fixture
def data():
    v = Vectors([[1, 2, 3], [1, 3, 1]])

    assert 2 == v.nb_vector
    assert 3 == v.nb_variable
    assert [1, 2] == v.get_identifiers()
    assert v
    return v

@pytest.fixture
def myi(data):
    return interface(data, "data/vectors.vec", Vectors)

@pytest.fixture
def vec10():
    return Vectors(get_shared_data("chene_sessile.vec"))



def test_identifiers():
    v = Vectors([[1, 2, 3], [4, 5, 6], [7, 8, 9]], Identifiers=[1, 2, 4])
    assert v.get_identifiers() == [1, 2, 4]

def test_types():
    v2 = Vectors([[1, 2.0, 3.0], [1, 5.0, 1.0]], Identifiers=[1, 2])

def test_constructor_from_file(myi):
    myi.constructor_from_file()

def test_constructor_from_file_failure(myi):
    myi.constructor_from_file_failure()

def test_constructor_one_variable():
    v = Vectors([1, 2, 3])
    assert v.nb_variable == 3

def test_constructor_identifiers_failure():
    try:
        # should be Vectors([[1,2,3], Identifiers=[1])
        v = Vectors([1, 2, 3], Identifiers=[1, 2])
        assert False
    except ValueError:
        assert True

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_vectors_pylist():
    """test vector constructor from list"""
    v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7]])
    assert v
    v = Vectors([[1.2, 0.34], [1.2, 0.34]])
    assert v
    v = Vectors([[1]])
    assert v
    v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1, 2, 3, 4]])
    assert v

    try:
        v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1, 2, 3]])
        assert False
    except:
        assert True

    try:
        v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7], [1.2, 2, 3]])
        assert False
    except:
        assert True

def test_len(data):
    v = data
    assert len(v) == 2
    assert len(v) == v.nb_vector

def test_plot(myi):
    myi.plot()

def _test_save(myi):
    myi.save(Format="Data")

def test_plot_write(myi):
    myi.plot_write()

def test_file_ascii_write(myi):
    myi.file_ascii_write()

def test_file_ascii_data_write(myi):
    myi.file_ascii_data_write()

def test_spreadsheet_write(myi):
    myi.spreadsheet_write()


def test_extract(data):
    """run and test the extract methods"""
    m = data
    m.extract(1)


def test_vectors_container():
    """vector container : len"""
    v = Vectors([[0, 1, 2, 3], [4, 5, 6, 7]])
    assert len(v) == 2

    for i in v:
        assert len(i) == 4

    assert v[0] == list(range(4))
    assert v[1][1] == 5

def test_variance_analysis(vec10):
    """test VarianceAnalysis"""
    # todo: finalise and make method variance_analysis robust.
    va = VarianceAnalysis(vec10, 1, 4, "O")
    assert vec10.variance_analysis(1, 4, 1, "whatever", variance_type["O"]) == str(
        va
    )

    try:
        va = VarianceAnalysis(vec10, 1, 4, "DUMMY")
        assert False
    except:
        assert True

def test_contingency_table(vec10):
    """test contingency table"""
    ct = ContingencyTable(vec10, 1, 4)
    assert ct and str(ct)

    ct2 = vec10.contingency_table(1, 4, "what", OutputFormat.ASCII)
    assert ct == ct2

def test_rank_computation(vec10):
    ComputeRankCorrelation(vec10, Type="Kendall", FileName="test")
    # ComputeRankCorrelation(vec10, Type="Spearman", FileName="test")

