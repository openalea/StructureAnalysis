""" Test renewal data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id:  $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.renewal import Renewal
from openalea.sequence_analysis.time_events import TimeEvents, NbEventSelect

from openalea.stat_tool.data_transform import * 
from openalea.sequence_analysis.data_transform import TimeScaling
from openalea.stat_tool.cluster import Cluster 
from openalea.stat_tool.cluster import Transcode, Cluster 

from tools import interface
from tools import runTestClass


class Test(interface):
    """to be done

    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "data/time_events.dat",
                           Renewal)
        
    def build_data(self):
        """todo: check identifier output. should be a list """
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        return TimeEvents('data/time_events.dat')
   
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
        #self.display_versus_str()
        
    def test_len(self):
        seq = self.data
        assert seq.nb_element == 42
        assert seq.nb_class == 8
        pass

    def _test_plot(self):        
        self.plot()
    
    def _test_save(self):
        self.save(skip_reading=True)
        self.save()
                    
    def test_plot_write(self):
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
    
    def _test_get_hmixture(self):
        data = self.data
        histo = data.get_hmixture()

    def test_nb_event_select(self):
        time = self.data
        res = time.nb_event_select(1, 4)
        assert res.nb_element == 7
        assert res.nb_class == 3

        res2 = NbEventSelect(time, 1, 4)
        assert str(res) == str(res2)

    def test_time_scaling(self):
        aml = self.data.time_scaling(2)
        mod = TimeScaling(self.data, 2)
        assert str(aml) == str(mod)
        
        
    def test_merge(self):
        time1 = self.data
        time2 = self.data
        assert str(Merge(time1,time2)) == str(time1.merge([time2]))

    def test_time_select(self):
        #max value must be greater than the offset.
        data= self.data
        data.time_select(3,35)



if __name__ == "__main__":
    runTestClass(Test())