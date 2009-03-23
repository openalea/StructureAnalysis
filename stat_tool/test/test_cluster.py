""" Cluster tests"""
__revision__ = "$Id: i$"


#from openalea.stat_tool import histogram, vectors, comparison, data_transform

from openalea.stat_tool.histogram import Histogram
from openalea.stat_tool.vectors import Vectors, VectorDistance
from openalea.stat_tool.comparison import Compare
from openalea.stat_tool.data_transform import SelectVariable
from openalea.stat_tool.cluster import Transcode, Clustering, ToDistanceMatrix, Cluster

class Test:
    """Test class to test cluster function and classes"""
    
    def __init__(self):
        pass

    def test_cluster(self):
        """test cluster on histograms"""

        fagus = Histogram("fagus1.his")

        histo2 = Cluster(fagus, "Step", 2)
        histo3 = Cluster(fagus, "Information", 0.8)
        histo4 = Cluster(fagus, "Limit", [2, 4, 6, 8, 10])

        assert histo2
        assert histo3
        assert histo4

    def test_transcode(self):
        """test transcode on histograms"""

        fagus = Histogram("fagus1.his")

        histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
        
        assert histo5

    
    def test_transcode_err(self):
        """test transcode errors on histograms"""

        fagus = Histogram("fagus1.his")

        try:
            _histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, ])
            assert False
        except:
            assert True

    def test_clustering(self):        
        """test clustering on vectors/matrices"""
        

        vec10 = Vectors("chene_sessile.vec")
        vec15 = SelectVariable(vec10, [1, 3, 6], Mode="Reject")

        assert vec15

        matrix10 = Compare(vec15, VectorDistance("N", "N", "N"))

        c1 = Clustering(matrix10, "Partition", 3, Prototypes=[1, 3, 12])
        Clustering(matrix10, "Hierarchy", Algorithm="Agglomerative")
        Clustering(matrix10, "Hierarchy", Algorithm="Divisive")

        assert c1
        assert ToDistanceMatrix(c1)





