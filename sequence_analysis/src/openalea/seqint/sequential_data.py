"""Stores generic sequential data type."""
from openalea.sequence_analysis import Sequences
from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences


class SequentialData(object):

    def __init__(self, input_sequence):
        self._input_sequence = input_sequence
        if input_sequence is not None:
            self._input_sequence = _MarkovianSequences(Sequences(self._input_sequence))
        self._restored_data = None

    @property
    def input_sequence(self):
        """input sequence getter."""
        return self._input_sequence

    @property
    def restored_data(self):
        """restored data getter."""
        return self._restored_data

    @restored_data.setter
    def restored_data(self, restored_data):
        """restored data setter."""
        self._restored_data = restored_data
