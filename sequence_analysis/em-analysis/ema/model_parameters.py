"""Stores model parameters."""


class ModelParameters(object):
    """Stores model parameters."""

    def __init__(self, model, initial_probabilities, transition_probabilities,
                 output_processes, occupancy_distributions):
        """Model Parameters constructor.

        :type model: Model
        :type initial_probabilities: InitialProbabilities
        :type transition_probabilities: TransitionProbabilities
        :type occupancy_distributions: list of OccupancyDistribution
        :type output_processes: list of OutputProcess

        :Example:
        >>> ModelParameters(model, initial_probabilities,
        transition_probabilities,
        [occupancy_distribution1, occupancy_distribution2],
        [output_process1, output_process2])
        """
        self._model = model
        self.initial_probabilities = initial_probabilities
        self.transition_probabilities = transition_probabilities
        self.occupancy_distributions = occupancy_distributions
        self.output_processes = output_processes

    def __str__(self):
        """Model parameters method for printing."""
        model_parameters_string = 'MODEL_PARAMETERS\n'
        model_parameters_string += str(self.initial_probabilities) + '\n'
        model_parameters_string += str(self.transition_probabilities) + '\n'
        for occupancy_distribution in self.occupancy_distributions:
            model_parameters_string += str(occupancy_distribution) + '\n'
        for output_process in self.output_processes:
            model_parameters_string += str(output_process) + '\n'
        return model_parameters_string

    def swap_hidden_state_order(self, k, i, j):
        assert k <= len(self.output_processes)-1
        assert i != j
        assert i <= self._model._k
        assert j <= self._model._k

        tmp = self.initial_probabilities._initial_probabilities[i]
        self.initial_probabilities._initial_probabilities[i] = self.initial_probabilities._initial_probabilities[j]
        self.initial_probabilities._initial_probabilities[j] = tmp

        tmp = self.transition_probabilities._transition_probabilities[i, :].copy()
        self.transition_probabilities._transition_probabilities[i, :] = self.transition_probabilities._transition_probabilities[j, :].copy()
        self.transition_probabilities._transition_probabilities[j, :] = tmp.copy()

        tmp = self.transition_probabilities._transition_probabilities[:, i].copy()
        self.transition_probabilities._transition_probabilities[:, i] = self.transition_probabilities._transition_probabilities[:, j].copy()
        self.transition_probabilities._transition_probabilities[:, j] = tmp.copy()

        tmp = self.occupancy_distributions[i]._state_number
        self.occupancy_distributions[i]._state_number = self.occupancy_distributions[j]._state_number
        self.occupancy_distributions[j]._state_number = tmp

        tmp = self.occupancy_distributions[i]
        self.occupancy_distributions[i] = self.occupancy_distributions[j]
        self.occupancy_distributions[j] = tmp

        tmp = self.output_processes[k]._observation_distributions[i]._observation_distribution_number
        self.output_processes[k]._observation_distributions[i]._observation_distribution_number = self.output_processes[k]._observation_distributions[j]._observation_distribution_number
        self.output_processes[k]._observation_distributions[j]._observation_distribution_number = tmp

        tmp = self.output_processes[k]._observation_distributions[i]
        self.output_processes[k]._observation_distributions[i] = self.output_processes[k]._observation_distributions[j]
        self.output_processes[k]._observation_distributions[j] = tmp

        self._model.update_hsmc_file()
