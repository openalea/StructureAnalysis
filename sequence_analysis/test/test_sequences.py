"""Distribution tests"""
__revision__ = "$Id: test_distribution.py 6219 2009-04-08 14:11:08Z cokelaer $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.sequences import Sequences

from tools import interface

class Test(interface):
    """a simple unittest class

    Integration Test 
    ================
    
    * 'ok' means works and testedPerform test on 
        
    ========================    ==================================
    ** from the interface**
    ascii_write                 ok
    display                     ok
    extract_data                nothing to be done
    file_ascii_write            ok
    plot                        ok                       
    save                        ok
    plot_print                  ok
    simulate                    ok
    plot_write                  ok
    spreadsheet_write           ok

    **others**
    old_plot                    ok   
    str                         ok
    len                         not relevant
    ========================    ==================================
 
    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           "sequences1.seq",
                           Sequences)
    
    def build_data(self):
        data = Sequences('sequences1.seq')
        return data
    
    def test_empty(self):
        self.empty()

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def test_print(self):
        #self.print_data()
        pass
        
    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()
        
    def test_len(self):
        len(self.data)
        
    
    def test_plot(self):        
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
        #self.simulate()
        pass
        
    def test_survival_ascii_write(self):
        pass
        #self.survival_ascii_write()
        
    def test_survival_spreadsheet_write(self):
        pass
        #self.survival_spreadsheet_write()
        
    def test_extract(self):
        pass

    def test_extract_data(self):
        pass 

   
