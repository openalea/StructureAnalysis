""" Regression tests"""
__version__ = "$Id: test_regression.py 9266 2010-07-11 13:02:34Z cokelaer $"


from openalea.stat_tool.regression import Regression
from openalea.stat_tool._stat_tool import _RegressionKernel as RegressionKernel
from openalea.stat_tool.vectors import Vectors
from openalea.stat_tool.compound import Compound
from openalea.stat_tool.distribution import Binomial

from tools import interface
from tools import runTestClass

class TestRegression(interface):
    """a simple unittest class"""

    def __init__(self):
        self.vector = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])
        interface.__init__(self,
                           self.build_data(),
                           None,
                           Regression)

    def build_data(self):

        #vector = Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])

        r1 = Regression(self.vector, "Linear", 1, 2)

        assert r1.nb_vector == 4
        return r1

    def test_build_bad_algorithm_failure(self):
        try:
            _r1 = Regression(self.vector, "Moving", 1, 2, 1,
                        Weighting=False)
            assert False
        except:
            assert True

    def test_get_residuals(self):
        for ii in range(0, self.data.nb_vector):
            assert self.data.get_residual(ii) == 0

    def test_print(self):
        self.print_data()

    def test_display(self):
        self.display()
        self.display_versus_ascii_write()
        self.display_versus_str()

    def test_len(self):
        """not implemented; irrelevant?"""
        assert self.data.nb_vector == 4

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
        pass

    def test_extract_data(self):
        pass

    def test_linear_regression(self):

        r1 = self.data

        #compare with the direct usage of linear regression
        r = self.vector.linear_regression(1, 2)

        assert r
        assert r1
        assert str(r) == str(r1)

    def test_moving_average(self):

        r1 = Regression(self.vector, "MovingAverage" , 1, 2, [1, ])
        r = self.vector.moving_average_regression_values(1, 2, [1, ], 'a')
        assert r
        assert r1
        assert str(r)==str(r1)

    def test_moving_average_failure(self):

        try:
            Regression(self.vector, "MovingAverage", 1, 2,  [1, ],
                       Algorithm="badAlgorithmName"
                       )
            assert False
        except:
            assert True

    def _test_moving_average_and_compound(self):
        """test to be implemented"""
        compound = Compound(Binomial(1, 10, 0.5), Binomial(1, 5, 0.4))
        Regression(self.vector, "MovingAverage", 1, 2, compound)

    def test_nearest_neighbours(self):

        r1 = Regression(self.vector, "NearestNeighbors", 1, 2, 1,
                        Weighting=False)
        r = self.vector.nearest_neighbours_regression(1, 2, 1., False)
        assert r
        assert r1
        assert str(r) == str(r1)

    def test_badtype(self):
        try:
            Regression(self.vector, "N", 1, 2, [1, ])
            assert False
        except TypeError:
            assert True


class _TestRegressionKernel():

    def __init__(self):
        self.data = RegressionKernel(4, 0, 10)

    def test_max_value(self):
        assert self.data.max_value == 10

    def test_min_value(self):
        assert self.data.min_value == 0

    def test_get_ident(self):
        assert self.data.ident == 4

    def others_to_be_done(self):
        #there are other methods that need to be tested with an
        #appropriate examples:
        #get_point
        #get_regression_df
        #get_residual_df
        #get_nb_parameter
        #get_parameter
        pass



if __name__ == "__main__":
    runTestClass(TestRegression())
