"""Hidden variable order markov tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id: test_hidden_variable_order_markov.py 9401 2010-08-10 12:24:59Z cokelaer $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.hidden_variable_order_markov import HiddenVariableOrderMarkov
from openalea.sequence_analysis.simulate import Simulate

from openalea.stat_tool.data_transform import *
from openalea.stat_tool.cluster import Cluster
from openalea.stat_tool.cluster import Transcode, Cluster

from tools import interface
from tools import runTestClass

from openalea.sequence_analysis import get_shared_data
\

class Test(interface):
    """a simple unittest class

    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           get_shared_data("dupreziana21.hc"),
                           HiddenVariableOrderMarkov)

    def build_data(self):
        """todo: check identifier output. should be a list """
        # build a list of 2 sequences with a variable that should be identical
        # to sequences1.seq
        hvom =  HiddenVariableOrderMarkov(get_shared_data('dupreziana21.hc'))

        return hvom

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

    def _test_simulate(self):
        sm = self.data
        sm.simulation_nb_elements(1, 10000, True)
        Simulate(sm,1, 10000, True)
        pass

    def _test_merge (self):
        s1 = self.data
        s2 = self.data

        sall = s1.merge([s2])

        #assert sall.nb_sequence == 4
        #assert sall.nb_variable == 2


    def test_thresholding(self):
        self.data.thresholding(1)

    def test_extract(self):
        #self.data.extract(1,0,0)
        pass

    def test_extract_data(self):
        pass
        #self.data.extract_data()



if __name__ == "__main__":
    runTestClass(Test())
