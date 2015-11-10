"""Semi markov tests

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr
"""
__revision__ = "$Id: test_semi_markov.py 8204 2010-02-19 10:27:45Z cokelaer $"


#from openalea.stat_tool import _stat_tool
#from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.nonhomogeneous_markov import NonhomogeneousMarkov
from openalea.sequence_analysis import get_shared_data
#from openalea.sequence_analysis.simulate import Simulate
#from openalea.sequence_analysis.sequences import Sequences
#from openalea.stat_tool.data_transform import *
#from openalea.stat_tool.cluster import Cluster
#from openalea.stat_tool.cluster import Transcode, Cluster

from tools import interface
from tools import runTestClass


def NonhomogeneousMarkovData():
    seq =  Sequences(get_shared_data('vanille_m.seq'))
    mc_m = Estimate(seq_m, "NONHOMOGENEOUS_MARKOV", "MONOMOLECULAR", "VOID")
    return mc_m

class Test(interface):
    """a simple unittest class for nonhomogeneous data

    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           get_shared_data("test_nonhomogeneous.dat"),
                           NonhomogeneousMarkov)

    def build_data(self):
        sm =  NonhomogeneousMarkov(get_shared_data('test_nonhomogeneous.dat'))
        return sm

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
        self.display_versus_str()

    def test_len(self):
        pass

    def test_plot(self):
        self.plot()

    def _test_save(self):
        self.save(skip_reading=True)

    def test_plot_write(self):
        self.plot_write()

    def _test_file_ascii_write(self):
        self.file_ascii_write()

    def _test_spreadsheet_write(self):
        self.spreadsheet_write()

    def _test_simulate(self):
        pass

    def test_extract(self):
        pass
        #self.data.extract(0,1)

    def test_extract_data(self):
        pass


#if __name__ == "__main__":
#    runTestClass(Test())
