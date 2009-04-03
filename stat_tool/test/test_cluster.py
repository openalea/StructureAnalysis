""" Cluster tests"""
__revision__ = "$Id: $"


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
        
        assert str(fagus.cluster_step(2))==str(histo2)
        assert str(fagus.cluster_information(0.8))==str(histo3)
        assert str(fagus.cluster_limit([2,4,6,8,10]))==str(histo4)



    def test_transcode(self):
        """test transcode on histograms"""

        fagus = Histogram("fagus1.his")

        histo5 = Transcode(fagus, [1, 2, 2, 3, 3, 4, 4, 5])
        
        assert str(histo5)==str(fagus.transcode([1, 2, 2, 3, 3, 4, 4, 5]))

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
        c2==matrix10.hierarchical_clustering(0, 2, "test", "test")
        c3==matrix10.hierarchical_clustering(0, 1, "test", "test")
        c4==matrix10.hierarchical_clustering(0, 0, "test", "test")
        
        # 1 for initialisation and 1 for divisive
        str(c1)==str(matrix10.partitioning_prototype(3, [1, 3, 12], 1, 1))
        str(c1_bis)==str(matrix10.partitioning_prototype(3, [1, 3, 12], 1, 2))

        #todo partioning_clusters
        

if __name__=="__main__":
    # perform all the test in the class Test (unit tests)
    test = Test()
    for method in dir(test):
        if method.startswith('_'):
            continue
        if callable(getattr(test, method)):
            getattr(test, method)()
        else:
            print 'skipping'
    # and functional tests.    




