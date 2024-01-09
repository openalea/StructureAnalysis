""" Test tops data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id: test_tops.py 9885 2010-11-06 18:19:34Z cokelaer $"

from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis import *

from tools import interface
from tools import runTestClass


def TopsData():
    """Returns simulated top"""
    param1 = TopParameters(get_shared_data("test_param1.p"), MaxPosition=20)
    tops = Simulate(param1, 100, 30)
    return tops


class Test(interface):
    """a simple unittest class


    """
    def __init__(self):
        interface.__init__(self,
                           self.build_data(),
                           get_shared_data("test_tops1.dat"),
                           Tops)

    def build_data(self):

        return Tops(get_shared_data('test_tops1.dat'))

    def _test_empty(self):
        self.empty()

    def test_constructor_from_file(self):
        self.constructor_from_file()

    def test_constructor_from_file_failure(self):
        self.constructor_from_file_failure()

    def _test_constructor_sequence(self):
        pass

    def _test_constructor_arrayn(self):
        pass

    def _test_constructor_array1(self):
        #print '--------------'
        #top = Tops([1,2,3,4,5], Identifiers=[1])
        print '------------'
        #print type(top)
        print '---------'
        print top
        #print self.data
        #top = Tops('data/tops1.dat')
        #print top

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_len(self):
        seq = self.data
        assert len(seq) == 2
        assert len(seq) == seq.nb_sequence

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

    def test_extract(self):
        """test to be done"""
        pass

    def test_extract_data(self):
        """test to be done"""
        pass

    def test_reverse(self):
        """test to be done"""
        pass

    def test_remove_apical_internodes(self):
        """test to be done"""
        pass


"""
t.display                              t.moving_average_from_distribution
t.estimation                           t.nb_sequence
t.estimation2                          t.nb_variable
t.extract                              t.old_plot
t.extract_length                       t.partial_autocorrelation_computation
t.extract_renewal_data                 t.plot
t.extract_sequence_length              t.plot_print
t.extract_time_events                  t.plot_write
t.extract_value                        t.pointwise_average
t.extract_vectors                      t.recurrence_time_sequences
t.file_ascii_data_write                t.remove_index_parameter
t.file_ascii_write                     t.remove_run
t.get_axillary_nb_internode            t.reverse
t.get_identifiers                      t.round
t.get_index_parameter_type             t.save
t.get_length                           t.scaling
t.get_max_value                        t.segmentation_array
t.get_min_value                        t.segmentation_change_point
t.alignment                            t.get_nb_internode                     t.segmentation_extract
t.alignment_vector_distance            t.get_plotable                         t.segmentation_model
t.ascii_data_write                     t.get_top_parameters                   t.segmentation_vector_distance
t.ascii_write                          t.get_type                             t.select_individual
t.build_nb_internode_histogram         t.index_parameter_extract              t.select_variable
t.build_vectors                        t.index_parameter_select               t.shift
t.cluster_limit                        t.length_select                        t.sojourn_time_sequences
t.cluster_step                         t.markovian_sequences                  t.spreadsheet_write
t.correlation_computation              t.max_length                           t.transcode
t.cross                                t.max_position                         t.transform_position
t.cumul_length                         t.merge                                t.value_select
t.cumulate                             t.merge_variable
t.difference                           t.moving_average
"""


if __name__ == "__main__":
    runTestClass(Test())
