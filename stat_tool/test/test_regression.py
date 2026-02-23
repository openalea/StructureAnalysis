"""Regression tests"""

try:
    from .tools import interface
    from .tools import robust_path as get_shared_data
except ImportError:
    from tools import interface
    from tools import robust_path as get_shared_data

from openalea.stat_tool._stat_tool import _RegressionKernel as RegressionKernel

from openalea.stat_tool.compound import Compound
from openalea.stat_tool.distribution import Binomial
from openalea.stat_tool.regression import Regression
from openalea.stat_tool.vectors import Vectors

import pytest

@pytest.fixture
def vector():
    return Vectors([[0, 0], [1, 1], [2, 2], [3, 3]])

@pytest.fixture
def data(vector):
    r1 = Regression(vector, "Linear", 1, 2)

    assert r1.nb_vector == 4
    return r1

@pytest.fixture
def myi(data):
    return interface(data, None, Regression)


def test_build_bad_algorithm_failure(vector):
    try:
        _r1 = Regression(vector, "Moving", 1, 2, 1, Weighting=False)
        assert False
    except:
        assert True

def test_get_residuals(data):
    for ii in range(0, data.nb_vector):
        assert data.get_residual(ii) == 0

def test_print(myi):
    myi.print_data()

def test_display(myi):
    myi.display()
    myi.display_versus_ascii_write()
    myi.display_versus_str()

def test_len(data):
    """not implemented; irrelevant?"""
    assert data.nb_vector == 4

def test_plot(myi):
    myi.plot()

def test_save(myi):
    myi.save(skip_reading=True)

def test_plot_write(myi):
    myi.plot_write()

def test_file_ascii_write(myi):
    myi.file_ascii_write()

def test_spreadsheet_write(myi):
    myi.spreadsheet_write()


def test_linear_regression(vector, data):
    r1 = data

    # compare with the direct usage of linear regression
    r = vector.linear_regression(1, 2)

    assert r
    assert r1
    #assert str(r) == str(r1)

def test_moving_average(vector):
    r1 = Regression(
        vector,
        "MovingAverage",
        1,
        2,
        [
            1,
        ],
    )
    from openalea.stat_tool.enums import algo_map

    r = vector.moving_average_regression_values(
        1,
        2,
        [
            1,
        ],
        algo_map["Averaging"],
    )
    assert r
    assert r1
    #assert str(r) == str(r1)

def test_moving_average_failure(vector):
    try:
        Regression(
            vector,
            "MovingAverage",
            1,
            2,
            [
                1,
            ],
            Algorithm="badAlgorithmName",
        )
        assert False
    except:
        assert True

def _test_moving_average_and_compound(vector):
    """test to be implemented"""
    compound = Compound(Binomial(1, 10, 0.5), Binomial(1, 5, 0.4))
    Regression(vector, "MovingAverage", 1, 2, compound)

def test_nearest_neighbours(vector):
    r1 = Regression(vector, "NearestNeighbors", 1, 2, 1, Weighting=False)
    r = vector.nearest_neighbours_regression(1, 2, 1.0, False)
    assert r
    assert r1
    #assert str(r) == str(r1)

def test_badtype(vector):
    try:
        Regression(
            vector,
            "N",
            1,
            2,
            [
                1,
            ],
        )
        assert False
    except TypeError:
        assert True


class _TestRegressionKernel:
    def __init__(self):
        self.data = RegressionKernel(4, 0, 10)

    def test_max_value(self):
        assert self.data.max_value == 10

    def test_min_value(self):
        assert self.data.min_value == 0

    def test_get_ident(self):
        assert self.data.ident == 4

    def others_to_be_done(self):
        # there are other methods that need to be tested with an
        # appropriate examples:
        # get_point
        # get_regression_df
        # get_residual_df
        # get_nb_parameter
        # get_parameter
        pass

