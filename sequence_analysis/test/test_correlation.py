""" Tests on ComputeAutoCorrelation, ComputeParialAutoCorrelation, 
ComputewhiteNoiseCorrelation 

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"


from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.data_transform import MovingAverage, VariableScaling
from openalea.sequence_analysis.correlation import ComputeCorrelation
from openalea.sequence_analysis.correlation import type_dict
from openalea.sequence_analysis.correlation import ComputeWhiteNoiseCorrelation
from openalea.sequence_analysis.correlation import ComputePartialAutoCorrelation

from openalea.stat_tool.distribution import Distribution


class Data():
    
    def __init__(self):
        
        self.sequence = self.create_sequence_data()
        self.type_map = type_dict
        
    def create_sequence_data(self):
        
        seq66 = Sequences("data/laricio_date66.seq")
        seq69 = MovingAverage(VariableScaling(seq66, 3, 100), 
                          Distribution("B", 0, 6, 0.5), BeginEnd=True, 
                          Output="Residual")
        return seq69
    

class test_ComputeCorrelation(Data):
        
    def __init__(self):
        Data.__init__(self)
        self.variable = 2
    def compute_correlation_type(self, variable, type, MaxLag=10, 
                                 Normalization="Exact"):
        seq = self.sequence
        cf = ComputeCorrelation(seq, variable, 
                               Type=type, MaxLag=MaxLag, Normalization=Normalization)
        assert cf.type == self.type_map[type]
        return cf
        
    def test_correlation_no_optional_arguments(self):
        seq = self.sequence
        cf = ComputeCorrelation(seq, self.variable)
        
    def _test_spearman(self):
        self.compute_correlation_type(self.variable, "Spearman", 10)
            
    def test_pearson(self):
        # used by test_ComputewhiteNoiseCorrelation
        return self.compute_correlation_type(self.variable, "Pearson", 10)
         
    def test_norm1(self):
        self.compute_correlation_type(self.variable, "Pearson", 10, "Approximated")
    
    def test_norm2(self):
        self.compute_correlation_type(self.variable, "Pearson", 10, "Exact")
    
    def test_norm3(self):
        try:
            self.compute_correlation_type(self.variable, "Pearson", 10, "Typolabel")
            assert False
        except Exception:
            assert True
    
            

class test_ComputeWhiteNoiseCorrelation(test_ComputeCorrelation):
    
    def __init__(self):
        test_ComputeCorrelation.__init__(self)
        self.correlation = self.test_pearson()
        
    def test_filter(self):
        data = self.correlation        
        ComputeWhiteNoiseCorrelation(data, [1, 1, 1])
        
    def test_order(self):
        data = self.correlation
        ComputeWhiteNoiseCorrelation(data, 1)
     
    def test_distribution(self): 
        data = self.correlation
        ComputeWhiteNoiseCorrelation(data , Distribution("Binomial", 0,4,0.5))
     
class test_ComputePartialAutoCorrelation(Data):
    
    def __init__(self):
        Data.__init__(self)
        
    def test_compute_partial_auto_correlation(self):
        ComputePartialAutoCorrelation(self.sequence, 2, MaxLag=5)
     
