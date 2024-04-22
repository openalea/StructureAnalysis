"""Stores occupancy distribution."""
import warnings


DEFAULT_DISTRIBUTION_NAME = 'NEGATIVE_BINOMIAL'
DEFAULT_STATE_TYPE = 'RECURRENT'
DEFAULT_LOWER_BOUND = 1
DEFAULT_UPPER_BOUND = None
DEFAULT_PARAMETER = 1
DEFAULT_PROBABILITY = 0.5


class OccupancyDistribution(object):
    """Stores occupancy distribution."""

    def __init__(self, model, state_number, state_type,
                 distribution_name=None, lower_bound=None, upper_bound=None,
                 parameter=None, probability=None):
        """Class constructor. Currently handle only discrete distributions.

        :param distribution_name: BINOMIAL, NEGATIVE_BINOMIAL, POISSON
        :type model: Model
        :type state_numer: int
        :type state_type: str in ['TRANSITORY', 'RECURRENT', 'ABSORBING']
        :type distribution_name: string in ['BINOMIAL', 'NEGATIVE_BINOMIAL',
            'POISSON']
        :type lower_bound: int
        :type upper_bound: int
        :type parameter: float
        :type probability: float

        :Example:
        >>> OccupancyDistribution(model, 0, 'RECURRENT',
            distribution_name='NEGATIVE_BINOMIAL', lower_bound=1,
            parameter=5, probability=0.3)

        >>> OccupancyDistribution(model, 0, 'ABSORBING')
        """
        call_error = True
        self._model = model
        self._state_number = state_number
        self._state_type = state_type
        self._distribution_name = distribution_name
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._parameter = parameter
        self._probability = probability

        if self._state_type == 'ABSORBING' and \
           all(item is None for item in
               [self._distribution_name, self._lower_bound, self._upper_bound,
                self._parameter, self._probability]):
            call_error = False
        elif self._distribution_name == 'BINOMIAL' and \
             self._parameter is None and all(item is not None for item in
             [self._lower_bound, self._upper_bound, self._probability]):
            call_error = False
        elif self._distribution_name == 'NEGATIVE_BINOMIAL' and \
             self._upper_bound is None and \
             all(item is not None for item in
             [self._lower_bound, self._parameter, self._probability]):
            call_error = False
        elif self._distribution_name == 'POISSON' and \
             all(item is not None for item in
             [self._lower_bound, self._parameter]) and \
             all(item is None for item in
             [self._upper_bound, self._probability]):
            call_error = False
        if call_error:
            raise AttributeError('Forbidden call to constructor.')

    def __str__(self):
        """Occupancy distribution method for printing."""
        return 'OCCUPANCY_DISTRIBUTION\n' \
               + 'state_number : ' + str(self._state_number) + '\n' \
               + 'state_type : ' + self._state_type + '\n' \
               + 'distribution_name : ' + str(self._distribution_name) + '\n' \
               + 'lower_bound : ' + str(self._lower_bound) + '\n' \
               + 'upper_bound : ' + str(self._upper_bound) + '\n' \
               + 'parameter : ' + str(self._parameter) + '\n' \
               + 'probability : ' + str(self._probability) + '\n' \


    @property
    def state_number(self):
        """state number getter."""
        return self._state_number

    @state_number.setter
    def state_number(self, state_number):
        """state number setter."""
        self._state_number = state_number
        self._model.update_hsmc_file()

    @state_number.deleter
    def state_number(self):
        """state number deleter."""
        del self._state_number

    @property
    def state_type(self):
        """state type getter."""
        return self._state_type

    @state_type.setter
    def state_type(self, state_type):
        """state type setter."""
        if state_type == 'ABSORBING':
            self._distribution_name = None
            self._lower_bound = None
            self._upper_bound = None
            self._parameter = None
            self._probability = None
            self._model.parameters.transition_probabilities._transition_probabilities[self._state_number, :] = \
                [0] * self._model._k
            self._model.parameters.transition_probabilities._transition_probabilities[
                self._state_number, self._state_number] = 1
        else:
            warnings.warn(
                'Occupancy distribution ' + str(self._state_number) + ' have been reset to default (GEOMETRIC p=0.5).',
                UserWarning)
            self._set_default_occupancy_parameters()
        self._state_type = state_type
        self._model.update_hsmc_file()

    @state_type.deleter
    def state_type(self):
        """state type deleter."""
        del self._state_type

    @property
    def distribution_name(self):
        """distribution name getter."""
        return self._distribution_name

    @distribution_name.setter
    def distribution_name(self, distribution_name):
        """distribution name setter."""
        if distribution_name == 'BINOMIAL':
            self._set_default_binomial_parameters()
        elif distribution_name == 'NEGATIVE_BINOMIAL':
            self._set_default_negative_binomial_parameters()
        elif distribution_name == 'POISSON':
            self._set_default_poisson_parameters()
        elif distribution_name == 'GEOMETRIC':
            self._set_default_binomial_parameters()
        else:
            warnings.warn('Unkown distribution name', UserWarning)
            return

        warnings.warn(
            'This state distribution has been changed. Think of'
            ' updating the parameters suited to the new distribution.',
            UserWarning)
        self._model.update_hsmc_file()

    @distribution_name.deleter
    def distribution_name(self):
        """distribution name deleter."""
        del self._distribution_name

    @property
    def lower_bound(self):
        """lower bound deleter."""
        return self._lower_bound

    @lower_bound.setter
    def lower_bound(self, lower_bound):
        """lower bound setter."""
        self._lower_bound = lower_bound
        self._model.update_hsmc_file()

    @lower_bound.deleter
    def lower_bound(self):
        """lower bound deleter."""
        del self._lower_bound

    @property
    def upper_bound(self):
        """upper bound deleter."""
        return self._upper_bound

    @upper_bound.setter
    def upper_bound(self, upper_bound):
        """upper bound setter."""
        self._upper_bound = upper_bound
        self._model.update_hsmc_file()

    @upper_bound.deleter
    def upper_bound(self):
        """upper bound deleter."""
        del self._upper_bound

    @property
    def parameter(self):
        """parameter getter."""
        return self._parameter

    @parameter.setter
    def parameter(self, parameter):
        """parameter setter."""
        self._parameter = parameter
        self._model.update_hsmc_file()

    @parameter.deleter
    def parameter(self):
        """parameter deleter."""
        del self._parameter

    @property
    def probability(self):
        """probability getter."""
        return self._probability

    @probability.setter
    def probability(self, probability):
        """probability setter."""
        self._probability = probability
        self._model.update_hsmc_file()

    @probability.deleter
    def probability(self):
        """probability deleter."""
        del self._probability

    def _set_default_occupancy_parameters(self):
        """reset occupancy distribution parameters to default."""
        self._distribution_name = DEFAULT_DISTRIBUTION_NAME
        self._state_type = DEFAULT_STATE_TYPE
        self._lower_bound = DEFAULT_LOWER_BOUND
        self._upper_bound = DEFAULT_UPPER_BOUND
        self._parameter = DEFAULT_PARAMETER
        self._probability = DEFAULT_PROBABILITY

    def _set_default_binomial_parameters(self):
        self._distribution_name = 'BINOMIAL'
        self._state_type = DEFAULT_STATE_TYPE
        self._lower_bound = DEFAULT_LOWER_BOUND
        self._upper_bound = 2
        self._parameter = 2
        self._probability = DEFAULT_PROBABILITY

    def _set_default_negative_binomial_parameters(self):
        self._set_default_occupancy_parameters()

    def _set_default_poisson_parameters(self):
        self._distribution_name = 'POISSON'
        self._state_type = DEFAULT_STATE_TYPE
        self._lower_bound = DEFAULT_LOWER_BOUND
        self._upper_bound = DEFAULT_UPPER_BOUND
        self._parameter = DEFAULT_PARAMETER
        self._probability = None
