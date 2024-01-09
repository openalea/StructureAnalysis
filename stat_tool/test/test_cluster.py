""" Cluster tests

:Author: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__version__ = "$Id: test_cluster.py 9380 2010-08-06 17:53:45Z cokelaer $"


from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.cluster import Transcode, Clustering, \
    ToDistanceMatrix, Cluster

from tools import runTestClass
from openalea.stat_tool import get_shared_data


class Test:
    """Test class to test cluster function and classes"""

    def __init__(self):
        self.data = self.build_data()

    def build_data(self):
        data = Histogram("data/fagus1.his")
        return data

    def test_cluster_histo(self):
        """test cluster on histograms"""
        fagus = self.data

        histo2 = Cluster(fagus, "Step", 2)
        histo3 = Cluster(fagus, "Information", 0.8)
        histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])

        assert histo2
        assert histo3
        assert histo4

        assert str(fagus.cluster_step(2)) == str(histo2)
        assert str(fagus.cluster_information(0.8)) == str(histo3)
        assert str(fagus.cluster_limit([2, 4, 6, 8, 10])) == str(histo4)

    def test_cluster_histo_failure(self):
        fagus = self.data

        try:
            _histo2 = Cluster(fagus, "Step", 2, 2)
            assert False
        except ValueError:
            assert True

    def test_cluster_vectors(self):
        v = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        assert str(Cluster(v, "Step", 1, 2)) == str(v.cluster_step(1, 2))
        assert str(Cluster(v, "Limit", 1, [2, 4, 6])) == \
            str(v.cluster_limit(1, [2, 4 ,6]))

    def _test_cluster_vectors1(self):
        v = Vectors([[1], [2], [3]])

        assert str(Cluster(v, "Step", 2)) == str(v.cluster_step(1, 2))
        assert str(Cluster(v, "Limit",  [2])) == \
            str(v.cluster_limit(1,[2]))

    def test_cluster_vectors_badtype(self):
        v = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        try:
            # should be step, cluster, information
            Cluster(v, "BadName", 2)
            assert False
        except KeyError:
            assert True

    def test_cluster_vectors_information_failure(self):
        v = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        try:
            # if v is vectors, information does not exist
            Cluster(v, "Information", 2)
            assert False
        except KeyError:
            assert True

    def test_cluster_discrete_sequence(self):
        """not yet implemented"""
        pass

    def test_transcode_histo(self):
        """test transcode on histograms"""
        fagus = self.data
        histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
        assert str(histo5)==str(fagus.transcode([1, 2, 2, 3, 3, 4, 4, 5]))

    def _test_transcode_vectors(self):
        vec = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        assert  str(vec.transcode(1, [1, 2, 3, 4]))==\
            str(Transcode(vec, 1, [1, 2, 3, 4]))

    def test_transcode_histo_err(self):
        fagus = self.data
        try:
            _histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, ])
            assert False
        except:
            assert True

    def test_clustering(self):
        vec10 = Vectors(get_shared_data("chene_sessile.vec"))
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")

        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))

        c1 = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12],
                        Algorithm="Divisive")
        c1_bis = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12],
                            Algorithm="Ordering")


        c2 = Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
        c3 = Clustering(matrix10, "Hierarchy", Algorithm="Divisive")
        c4 = Clustering(matrix10, "Hierarchy", Algorithm="Ordering")

        assert c1
        assert c2
        assert c3
        assert c4
        assert ToDistanceMatrix(c1)

        # first argument is the Algorithm
        #  * 0 for agglomerative
        #  * 1 for divisive
        # Second argument is the criterion
        #  * 2 for averaging

        #those 3 tests works on my laptop (TC, April 2009) but not on buildbot
        #assert c2 == matrix10.hierarchical_clustering(0, 2, "test", "test")
        #assert c3 == matrix10.hierarchical_clustering(1, 1, "test", "test")
        #assert c4 == matrix10.hierarchical_clustering(2, 0, "test", "test")

        # 1 for initialisation and 1 for divisive
        assert str(c1) == \
            str(matrix10.partitioning_prototype(3, [1, 3, 12], 1, 1))
        assert str(c1_bis) == \
            str(matrix10.partitioning_prototype(3, [1, 3, 12], 1, 2))



if __name__ == "__main__":
    runTestClass(Test())
