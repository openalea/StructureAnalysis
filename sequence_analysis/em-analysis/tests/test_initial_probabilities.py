"""Testing class InitialProbabilities."""
from ema import EyeMovementData
from ema import InitialProbabilities
from ema import Model
from ema import TEST_PATH
import os
import pytest


@pytest.fixture
def model():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def initial_probabilities():
    return InitialProbabilities(None, [0.1, 0.2, 0.5, 0.2, 0])


def test_initial_probabilities_constructor(initial_probabilities):
    assert((initial_probabilities._initial_probabilities == \
        [0.1, 0.2, 0.5, 0.2, 0]).all())


def test_initial_probabilities_getter(initial_probabilities):
    assert((initial_probabilities.initial_probabilities == \
        [0.1, 0.2, 0.5, 0.2, 0]).all())


def test_initial_probabilities_setter(model):
    model.parameters.initial_probabilities.initial_probabilities = \
        [0.1, 0.2, 0.3, 0.4, 0]
    assert((model.parameters.initial_probabilities.initial_probabilities == \
            [0.1, 0.2, 0.3, 0.4, 0]).all())
