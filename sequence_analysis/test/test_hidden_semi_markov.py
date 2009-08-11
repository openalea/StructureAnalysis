"""hidden semi markov data structure  tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

"""
__revision__ = "$Id:  $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.hidden_semi_markov import HiddenSemiMarkov
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.data_transform import Thresholding

from openalea.stat_tool.data_transform import * 
from openalea.stat_tool.cluster import Cluster 
from openalea.stat_tool.cluster import Transcode, Cluster 

from tools import interface
from tools import runTestClass


class Test(interface):
    """a simple unittest class
    
    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/hidden_semi_markov.dat",
                           HiddenSemiMarkov)
        
    def build_data(self):
        """todo: check identifier output. should be a list """
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        hsm = HiddenSemiMarkov('data/hidden_semi_markov.dat')
        return hsm
   
    def test_empty(self):
        pass
        #self.empty()

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
        
    def test_simulate(self):
        sm = self.data
        sm.simulation_nb_elements(1, 10000, True)
        Simulate(sm,1, 10000, True)
        pass
        
    def test_thresholding(self):
        a = self.data.thresholding(0.01)
        b = Thresholding(self.data, MinProbability=0.01)
        #assert str(a)==str(b) 

    def test_extract(self):
        self.data.extract(1,1,1)

    def test_extract_data(self):
        self.data.extract_data()
         

"""
hsm.nb_output_process
hsm.ascii_write                                            
hsm.divergence_computation          
hsm.simulation_histogram
hsm.extract_histogram               
hsm.simulation_markovian_sequences
hsm.file_ascii_write                
hsm.get_forward  
hsm.simulation_nb_sequences
hsm.get_plotable                    
hsm.get_semi_markov_data
hsm.state_sequence_computation
hsm.get_state_subtype               
hsm.nb_iterator                     
"""


if __name__ == "__main__":
    runTestClass(Test())
