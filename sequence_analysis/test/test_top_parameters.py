""" Test tops data structure

.. author:: Thomas Cokelaer, Thomas.Cokelaer@inria.fr

.. todo:: to be done
"""
__revision__ = "$Id: test_top_parameters.py 9885 2010-11-06 18:19:34Z cokelaer $"


from openalea.stat_tool import _stat_tool
from openalea.sequence_analysis import _sequence_analysis
from openalea.sequence_analysis.tops import Tops
from openalea.sequence_analysis.simulate import Simulate
from openalea.sequence_analysis.top_parameters import TopParameters
from openalea.sequence_analysis.data_transform import *
from openalea.sequence_analysis import get_shared_data



from tools import interface
from tools import runTestClass


def TopParametersData():
    """Returns simulated top"""
    param1 = TopParameters(get_shared_data("test_param1.p"), MaxPosition=20)
    return param1

class Test(interface):
    """a simple unittest class


    """
    def __init__(self):
        interface.__init__(self,
                           self.test_build_data(),
                           get_shared_data("test_param1.p"),
                           TopParameters)


    def test_build_data(self):
        return TopParameters(0.6, 0.6, 1.2)


    def _test_empty(self):
        #test tobedone
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

    def test_save(self):
        pass
        #self.save(skip_reading=True)

    def test_plot_write(self):
        self.plot_write()

    def _test_file_ascii_write(self):
        self.file_ascii_write()

    def test_spreadsheet_write(self):
        self.spreadsheet_write()

    def test_extract(self):
        """test to be done"""
        from openalea.stat_tool._stat_tool import _DiscreteParametricModel
        for i in range(1, self.data.max_position):
            assert type(self.data.extract(i)) == _DiscreteParametricModel

    def test_extract_data(self):
        """test to be done"""
        pass
    def test_axillary_probability(self):
        assert self.data.axillary_probability == 0.6

    def test_probability(self):
        assert self.data.probability == 0.6

    def test_rhythm_ratio(self):
        assert self.data.rhythm_ratio == 1.2

    def test_get_axillary_nb_internode(self):
        for i in range(1, self.data.max_position):
            self.data.get_axillary_nb_internode(i)

    def test_get_axillary_nb_internode_wrong_arg(self):

        try:
            self.data.get_axillary_nb_internode(0)
            assert False
        except:
            assert True

        try:
            self.data.get_axillary_nb_internode(self.data.max_position+1)
            assert False
        except:
            assert True

    def test_simulate(self):
        top1 = Simulate(self.data, 200, 30, NbAxillary=2)
        top2 = self.data.simulate(200, 30, 2)


"""
remain to be tested :
top.ascii_write
top.simulation_dists
top.get_tops

"""


if __name__ == "__main__":
    runTestClass(Test())
