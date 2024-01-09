"""Testing class TransitionProbabilities."""
from ema import EyeMovementData
from ema import Model
from ema import TransitionProbabilities
from ema import TEST_PATH
import numpy as np
import os
import pytest


@pytest.fixture
def model():
    return Model(EyeMovementData(), init_hsmc_file=os.path.join(TEST_PATH, 'input_model.hsmc'))


@pytest.fixture
def model_set_transition_probabilities(model):
    model.parameters.transition_probabilities.transition_probabilities = (
        np.matrix([[0, 0, 0, 0.5, 0.5], [0, 0, 0.5, 0.5, 0], [0.5, 0.5, 0, 0, 0],
                   [0, 0.5, 0.5, 0, 0], [0.5, 0, 0, 0.5, 0]]))
    return model


@pytest.fixture
def transition_probabilities():
    return TransitionProbabilities(
        None, np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))


def test_transition_probabilities_constructor(transition_probabilities):
    assert (transition_probabilities._transition_probabilities ==
            np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).all()


def test_transition_probabilities_getter(transition_probabilities):
    assert (transition_probabilities.transition_probabilities ==
            np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])).all()


def test_transition_probabilities_matrix_setter(model_set_transition_probabilities):
    assert (model_set_transition_probabilities.parameters.transition_probabilities.transition_probabilities ==
            np.matrix([[0, 0, 0, 0.5, 0.5], [0, 0, 0.5, 0.5, 0], [0.5, 0.5, 0, 0, 0],
                       [0, 0.5, 0.5, 0, 0], [0.5, 0, 0, 0.5, 0]])).all()


def test_transition_probabilities_line_setter(model_set_transition_probabilities):
    model_set_transition_probabilities.parameters.transition_probabilities.\
    transition_probabilities = ([0, 0, 1, 0, 0], 0)
    assert (model_set_transition_probabilities.parameters.transition_probabilities.transition_probabilities ==
            np.matrix([[0, 0, 1, 0, 0], [0, 0, 0.5, 0.5, 0], [0.5, 0.5, 0, 0, 0],
                       [0, 0.5, 0.5, 0, 0], [0.5, 0, 0, 0.5, 0]])).all()
