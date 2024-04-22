"""Stores observation distribution."""


class ObservationDistribution(object):
    """Stores occupancy distribution."""

    def __init__(self, model, observation_distribution_number, outputs):
        """Observation Distribution constructor.

        :type model: Model
        :type observation_distribution_number: int
        :type outputs: list of float

        :Example:
        >>> ObservationDistribution(model, 1, [0.2, 0.8])
        """
        self._model = model
        self._observation_distribution_number = observation_distribution_number
        self._outputs = outputs

    def __str__(self):
        """observation distribution method for printing."""
        return 'OBSERVATION_DISTRIBUTION\n'  \
               + 'observation_distribution_number : '  \
               + str(self._observation_distribution_number) + '\noutputs : '  \
               + str(self._outputs) + '\n'

    @property
    def observation_distribution_number(self):
        """observation distribution number getter."""
        return self._observation_distribution_number

    @observation_distribution_number.setter
    def observation_distribution_number(self, observation_distribution_number):
        """observation distribution number setter."""
        self._observation_distribution_number = observation_distribution_number
        self._model.update_hsmc_file()

    @observation_distribution_number.deleter
    def observation_distribution_number(self):
        """observation distribution number deleter."""
        del self._observation_distribution_number

    @property
    def outputs(self):
        """outputs getter."""
        return self._outputs

    @outputs.setter
    def outputs(self, value):
        """outputs setter.

        :type value: list of float or numpy.ndarray of float

        :Example:
        >>> outputs([0.2, 0.4, 0.6])
        """
        if all(isinstance(item, float) for item in value):
            for i in xrange(0, len(self._outputs)):
                self._outputs[i] = value[i]
        else:
            raise AttributeError("Forbidden usage.")
            return

        self._model.update_hsmc_file()

    @outputs.deleter
    def outputs(self):
        """outputs deleter."""
        del self._outputs
