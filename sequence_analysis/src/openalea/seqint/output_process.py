"""Stores output processes."""


class OutputProcess(object):
    """Stores output processes."""

    def __init__(self, model, output_process_number, output_process_type,
                 observation_distributions):
        """"Output Process constructor.

        :type model: Model
        :type output_process_number: int
        :type output_process_type: String in ['CATEGORICAL']
        :type observation_distributions: list of ObservationDistribution

        :Example:
        >>> OutputProcess(model, 0, 'CATEGORICAL', [observation_distribution1,
        observation_distribution2])
        """
        if output_process_type != 'CATEGORICAL':
            raise AttributeError('Forbidden call.')
            return
        self._model = model
        self._output_process_number = output_process_number
        self._output_process_type = output_process_type
        self._observation_distributions = observation_distributions

    def __str__(self):
        """Output Process method for printing."""
        output_process_string = 'OUTPUT_PROCESS\noutput_process_number : '
        output_process_string += str(self._output_process_number)
        output_process_string += '\noutput_process_type : '
        output_process_string += self._output_process_type
        output_process_string += '\n'
        for observation_distribution in self._observation_distributions:
            output_process_string += str(observation_distribution)
        return output_process_string

    @property
    def output_process_number(self):
        """output process number getter."""
        return self._output_process_number

    @output_process_number.setter
    def output_process_number(self, output_process_number):
        """output process number setter."""
        self._output_process_number = output_process_number
        self._model.update_hsmc_file()

    @output_process_number.deleter
    def output_process_number(self):
        """output process number deleter."""
        del self._output_process_number

    @property
    def output_process_type(self):
        """output process type getter."""
        return self._output_process_type

    @output_process_type.setter
    def output_process_type(self, output_process_type):
        """output process type setter."""
        if output_process_type is not 'CATEGORICAL':
            raise AttributeError('Forbidden call.')
            return
        self._output_process_type = output_process_type
        self._model.update_hsmc_file()

    @output_process_type.deleter
    def output_process_type(self):
        """output process type deleter."""
        del self._output_process_type

    @property
    def observation_distributions(self):
        """observation_distributions getter."""
        return self._observation_distributions
