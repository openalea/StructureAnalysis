"""Semi markov tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id$"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.variable_order_markov import VariableOrderMarkov
from openalea.sequence_analysis.sequences import Sequences
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.estimate import Estimate

from openalea.stat_tool.data_transform import * 
from openalea.stat_tool.cluster import Cluster 
from openalea.stat_tool.cluster import Transcode, Cluster 

from tools import interface
from tools import runTestClass


def VariableOrderMarkovData():
    sm =  VariableOrderMarkov('data/variable_order_markov.dat')
    ret = Simulate(sm, 1, 1000, True)
    return sm

    
class Test(interface):
    """a simple unittest class
    
    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/variable_order_markov.dat",
                           VariableOrderMarkov)
        
    def build_data(self):
        seq = Sequences("data/belren1.seq")
        vom = Estimate(seq, "VARIABLE_ORDER_MARKOV", "Ordinary", 
                        MaxOrder=4, GlobalInitialTransition=False)
                
        return vom
   
    def _test_empty(self):
        self.empty()

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_print(self):
        self.print_data()
        
    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()
        
    def test_len(self):
        seq = self.data
        pass

    def _test_plot(self):        
        self.plot()
    
    def test_save(self):
        self.save(skip_reading=True)
                    
    def test_plot_write(self):
        self.plot_write()
        
    def test_file_ascii_write(self):
        self.file_ascii_write()
        
    def test_spreadsheet_write(self):
        self.spreadsheet_write()
        
    def _test_simulate(self):
        sm = self.data
        #sm.simulation_nb_elements(1, 10000, True)
        Simulate(sm,1, 10000, True)
        pass
        
    def test_thresholding(self):
        self.data.thresholding(1)
        
    def test_extract(self):
        self.data.extract(1,0,0)

    def test_extract_data(self):
        self.data.extract_data()


if __name__ == "__main__":
    runTestClass(Test())