"""vectors tests"""

__version__ = "$Id$"

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data

from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.vectors import VectorDistance, Vectors
from openalea.stat_tool.cluster import Clustering

from pathlib import Path

import pytest

@pytest.fixture
def data():
    vec10 = Vectors(get_shared_data("chene_sessile.vec"))
    vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
    matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
    assert 138 == matrix10.nb_row
    assert 138 == matrix10.nb_column
    return matrix10

@pytest.fixture
def myi(data):
    # init expect the 4th argument to be provided.
    # vectors is therefore passed as dummy structure
    return interface(data, None, Vectors)

def test_get_column_identifier(data):
    assert data.get_column_identifier(0) == 1

def test_get_deletion_distance(data):
    assert data.get_deletion_distance(0,0) == -1.

def test_get_distance(data):
    assert data.get_distance(0,0) == 0.

def test_get_insertion_distance(data):
    assert data.get_inserion_distance(0,0) == -1.

def test_get_nb_match(data):
    assert data.get_nb_match(0,0) == -1

def test_get_length(data):
    assert data.get_length(0,0) == 0.0

def test_get_insertion(data):
    assert data.get_nb_insertion(0,0) == -1

def test_get_deletion(data):
    assert data.get_nb_deletion(0,0) == -1

def test_get_nb_substitution(data):
    assert data.get_nb_substitution(0,0) == -1

def test_get_nb_transposition(data):
    assert data.get_nb_transposition(0,0) == -1

def test_get_row_identifier(data):
    assert data.get_row_identifier(0) == 1

def test_get_substitution_distance(data):
    assert data.get_substitution_distance(0,0) == -1

def test_get_insertion_distance(data):
    assert data.get_insertion_distance(0,0) == -1

def test_get_shape(data):
    assert (data.nb_column, data.nb_row) == (138, 138)

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_plot(myi):
    myi.plot()

def test_plot_write(myi):
    myi.plot_write()

def test_select_individual(myi):
    assert myi.data.select_individual([1], False)
    
def test_select_individual(myi):
    assert myi.data.select_individual([1], False)
    
def test_wrong_partitioning_clusters(myi):
    try:
        myi.data.partitioning_clusters([0])
        myi.data.partitioning_clusters([])
        myi.data.partitioning_clusters([[1,2], [3,4]])
        assert False
    except:
        assert True

def test_partitioning_clusters(myi):
   clust1 = Clustering(myi.data, "Partition", 2)
   nb_clusters = clust1.get_nb_cluster()
   # partition
   part1 = [[i for i in range(1,myi.data.nb_row+1) if clust1.get_assignment(i) == c] for c in range(1,nb_clusters+1)]
   assert myi.data.partitioning_clusters(part1)

def test_wrong_partitioning_prototype(myi):
    """
    Testing errors in partitioning_prototype.
    Other features of partitioning_prototype are tested in test_cluster.py
    """
    try:
        assert myi.data.partitioning_prototype(0, [0], 0, 0)
        assert False
    except:
        assert True    
            

def test_wrong_symmetrize(myi):
    try:
        myi.data.symmetrize()
        assert False
    except:
        assert True    
    
def test_test_symmetry(myi):
    assert myi.data.test_symmetry()

def test_unnormalize(myi):
    try:
        myi.data.unnormalize()
        assert False
    except:
        assert True    
    
def test_spreadsheet_write(myi):
    myi.spreadsheet_write()


if __name__ == "__main__":
    def path():
        return Path(__file__).parent

    data_file = str(get_shared_data("chene_sessile.vec"))

    def data():
        vec10 = Vectors(get_shared_data(data_file))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")
        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))
        return matrix10
    
    
    
    myi = interface(data,  data_file, Vectors)
    myi.data = data()

    test_get_column_identifier(data())
    test_select_individual(myi)
    test_partitioning_clusters(myi)
    test_wrong_partitioning_clusters(myi)