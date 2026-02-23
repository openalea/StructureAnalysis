"""vectors tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
from .tools import interface, runTestClass

from openalea.stat_tool.vectors import VectorDistance

import pytest

@pytest.fixture
def data():
    v = VectorDistance("N", "O", "S")
    assert v.get_distance_type() == 0
    assert v
    return v


@pytest.fixture
def myi(data):
    return interface(data, "data/vector_distance.vd", VectorDistance)


def test_constructor_from_file(myi):
    myi.constructor_from_file()

def test_constructor_from_file_failure(myi):
    myi.constructor_from_file_failure()

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_len(data):
    v = data
    assert len(v) == 3
    assert len(v) == v.nb_variable

def test_save(myi):
    myi.save()

def test_file_ascii_write(myi):
    try:
        myi.file_ascii_write()
    except Exception as error:
        pass
        print(error)

def test_spreadsheet_write(myi):
    try:
        myi.spreadsheet_write()
    except Exception as error:
        pass
        print(error)


def test_vector_distance():
    """test vector distance constructors"""
    v = VectorDistance("NUMERIC", "ORDINAL", "SYMBOLIC")
    assert v and len(v) == 3

    v = VectorDistance(2.3, "N", 4, "O", 6, "S")
    assert v and len(v) == 3

    # wrapper should work. needs python interface ?
    # v = VectorDistance( (2.3, 'N'),  (4, 'O'), (6, 'S'))
    # assert v and len(v) == 3

    v = VectorDistance(2.3, "N", 4, "O", 6, "S", Distance="QUADRATIC")
    assert v and len(v) == 3

    v = VectorDistance("NUMERIC", "ORDINAL", "SYMBOLIC", Distance="QUADRATIC")
    assert v and len(v) == 3

    assert str(VectorDistance("N", "O", "S"))

def test_constructor_failure():
    try:
        v = VectorDistance("N", "DUMMY")
        assert False
    except:
        assert True

