"""Distribution tests"""
__revision__ = "$Id: test_distribution.py 6219 2009-04-08 14:11:08Z cokelaer $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.time_events import TimeEvents

from openalea.stat_tool.data_transform import * 
from openalea.stat_tool.cluster import Cluster 
from openalea.stat_tool.cluster import Transcode, Cluster 

from tools import interface

class Test(interface):
    """a simple unittest class

    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/time_events.dat",
                           TimeEvents)
        
    def build_data(self):
        """todo: check identifier output. should be a list """
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        return TimeEvents('data/time_events.dat')
   
    def test_empty(self):
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
        #self.display_versus_str()
        
    def test_len(self):
        seq = self.data
        assert seq.nb_element==42
        assert seq.nb_class==8
        pass

    def test_plot(self):        
        self.plot()
    
    def test_save(self):
        self.save(skip_reading=True)
        self.save()
                    
    def _test_plot_write(self):
        self.plot_write()
        
    def test_file_ascii_write(self):
        self.file_ascii_write()
        
    def test_spreadsheet_write(self):
        self.spreadsheet_write()
        
    def _test_simulate(self):
        #self.simulate()
        pass
        
        
    def test_extract(self):
        """todo"""
        pass 

    def test_extract_data(self):
        """todo"""
        pass 

    def test_get_htime(self):
        data = self.data
        histo = data.get_htime()


    def test_get_hnb_event(self):
        data = self.data
        histo = data.get_hnb_event(20)
    
    def test_get_hmixture(self):
        data = self.data
        histo = data.get_hmixture()
