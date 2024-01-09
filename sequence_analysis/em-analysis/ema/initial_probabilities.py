"""Stores initial probabilities."""
import numpy as np


class InitialProbabilities(object):
    """Stores oculometric data as pandas dataframe."""

    def __init__(self, model, initial_probabilities):
        """Initial Probabilities constructor.

        :type model: Model
        :type initial_probabilities: list of float or numpy.ndarray of float

        :Example:
        >>> InitialProbabilities(model, [0.1, 0.2, 0.5, 0.2])
        """
        self._model = model
        self._initial_probabilities = np.array(initial_probabilities)

    def __str__(self):
        """Initial Probabilities method for printing."""
        return 'INITIAL PROBABILITIES\n' \
               + str(self._initial_probabilities) \
               + '\n'

    @property
    def initial_probabilities(self):
        """initial probabilities getter."""
        return self._initial_probabilities

    @initial_probabilities.setter
    def initial_probabilities(self, value):
        """initial proabilities setter.

        If value is a list then change the entire list.
        If value is an int or float and index is specified then set the
        value specified at index index.

        :type value: int/float or list of int/float
        :type index: int

        :Example:
        >>> initial_probabilities([0.6, 0.4])
        """
        if any(isinstance(value, item) for item in [list, np.ndarray]):
            self._initial_probabilities = np.array(value)

        else:
            raise ValueError('Forbidden usage.')
            return

        self._model.update_hsmc_file()

    @initial_probabilities.deleter
    def initial_probabilities(self):
        """initial probabilities deleter."""
        del self._initial_probabilities
