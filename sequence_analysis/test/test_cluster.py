""" Cluster tests"""
__revision__ = "$Id: test_cluster.py 6325 2009-04-29 16:20:55Z cokelaer $"


from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.semi_markov import SemiMarkov
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.histogram import Histogram


class _Cluster:
    """Test class to test cluster function and classes"""
    
    def __init__(self):
        self.data = None    
    
    def create_data(self):
        raise NotImplemented
        
    def test_cluster_step(self):
        raise NotImplemented
    
    def test_cluster_limit(self):
        raise NotImplemented

class _HistoCase():
    def test_cluster_information(self):
            raise NotImplemented
    
class TestHistogram(_Cluster, _HistoCase):
    
    def __init__(self):
        _Cluster.__init__(self)
        self.data = self.create_data()
    
    def create_data(self):
        return Histogram('data/fagus1.his')
    
    def test_cluster_step(self):
        cluster1 = Cluster(self.data, "Step", 2)
        cluster2 = self.data.cluster_step(2)
        assert str(cluster1) == str(cluster2)
        
    def test_cluster_limit(self):
        cluster1 = Cluster(self.data, "Limit", [2, 4, 6, 8, 10])
        cluster2 = self.data.cluster_limit([2, 4, 6, 8, 10])
        assert str(cluster1) == str(cluster2)

    def test_cluster_information(self):
        cluster1 = Cluster(self.data, "Information", 0.8)
        cluster2 = self.data.cluster_information(0.8)
        assert str(cluster1) == str(cluster2)

#    def test_cluster_vectors(self):
        #v = Vectors([[1, 2, 3], [1, 3, 1], [4, 5, 6]])
        #assert str(Cluster(v, "Step", 1, 2)) == str(v.cluster_step(1, 2))  
        #assert str(Cluster(v, "Limit", 1, [2, 4, 6])) == \
            #str(v.cluster_limit(1, [2, 4 ,6]))        


