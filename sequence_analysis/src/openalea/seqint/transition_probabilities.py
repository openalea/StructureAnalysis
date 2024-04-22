"""Stores transition probabilities."""
import numpy as np
import pytest

class TransitionProbabilities(object):
    """Stores oculometric data as pandas dataframe."""

    def __init__(self, model, transition_probabilities):
        """Transition Probabilities constructor.

        :type model: Model
        :type transition_probabilities: numpy.matrixlib.defmatrix.matrix

        :Example:
        >>> TransitionProbabilities(model, np.matrix([[1, 2, 3], [4, 5, 6],
        [7, 8, 9]]))
        """
        self._model = model
        self._transition_probabilities = transition_probabilities

    def __str__(self):
        """Transition Probabilities method for printing."""
        return 'TRANSITION_PROBABILITIES\n' \
               + str(self._transition_probabilities) + '\n'

    @property
    def transition_probabilities(self):
        """transition probabilities getter."""
        return self._transition_probabilities

    @transition_probabilities.setter
    def transition_probabilities(self, value):
        """transition probabilities setter.

        If value is a numpy matrix then change the entire matrix
        If value is a list and line_index is specified
        then change it.

        :type value: tuple(list of int/float or numpy.ndarray of
        int/float, line_index) OR numpy.matrixlib.defmatrix.matrix of int/float

        :Example:
        >>> transition_probabilities = ([0.2, 0.2, 0.6], 3)

        >>> transition_probabilities =
                np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        """

        if isinstance(value, np.matrixlib.defmatrix.matrix):
            for i in xrange(0, value.shape[0]):
                assert(pytest.approx(np.sum(value[i, :])) == 1)
                if value[i, i] == 1:
                    self._model.parameters.occupancy_distributions[
                        i].state_type = 'ABSORBING'
                    print('STATE ' + str(i) +
                          ' has been set to absorbing.')
                    return
            self._transition_probabilities = value

        elif (isinstance(value, tuple) and any(isinstance(value[0], item) \
             for item in [list, np.ndarray])):
            prob, line_index = value
            if prob[line_index] == 1:
                self._model.parameters.occupancy_distributions[line_index].\
                    state_type = 'ABSORBING'
                print('STATE ' + str(line_index) +
                      ' has been set to absorbing.')
                return
            assert(pytest.approx(sum(prob)) == 1)
            self._transition_probabilities[line_index, :] = prob

        else:
            raise ValueError('Forbidden usage.')
            return

        self._model.update_hsmc_file()

    @transition_probabilities.deleter
    def transition_probabilities(self):
        """transition probabilities deleter."""
        del self._transition_probabilities
