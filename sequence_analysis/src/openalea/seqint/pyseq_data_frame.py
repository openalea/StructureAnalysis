"""
Represent a set of sequences as a pandas.DataFrame
"""
import math
from openalea.sequence_analysis import Sequences
from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences
import pandas as pd

class PySeqDataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        """
        Init a PySeqDataFrame.
        seq_index is the number of the column containing the index of sequences
        """
        index_as_ms = kwargs.pop('index_as_ms', False)
        seq_index_name = kwargs.pop('seq_index_name', False)
        super(PySeqDataFrame, self).__init__(*args, **kwargs)
        self.index_as_ms = index_as_ms if type(index_as_ms) == bool else None
        self.seq_index_name = seq_index_name if type(seq_index_name) == str else None
        self._input_sequence = None
        self._input_sequence_getter_optimizer_last = None
        self._input_sequence_getter_optimizer_bool = False
        self.restored_data = None

    def col_to_seq(self, col_names):
        """Extract cols from DataFrame as a sequence of sequences."""
        if all(isinstance(item, str) for item in col_names):
            for col_name in col_names:
                assert(col_name in self.columns)
        else:
            raise ValueError('All %s must be in %s' % (col_names, self.column_names))
        data = []
        is_first = self[self.seq_index_name].diff() > 0
        is_first[1] = True
         # End of first sequence
        seq_id = 1
        for i in self.index[is_first]:    
            if i < max(self.index[is_first]):
                seq_end_index = self.index[is_first][seq_id]
                data.append(self.loc[range(i, seq_end_index), col_names].values.tolist())
                seq_id += 1
        seq_end_index = max(self.index)+1
        data.append(self.loc[range(i, seq_end_index), col_names].values.tolist())
        return data

    def get_input_sequence(self, col_names):
        """_input_sequence getter using optimizers based on last request"""
        if (not self._input_sequence_getter_optimizer_bool) or (
                self._input_sequence_getter_optimizer_last != col_names):
            self._input_sequence_getter_optimizer_bool = True
            self._input_sequence = _MarkovianSequences(Sequences(self.col_to_seq(col_names)))
        else:
            self._input_sequence_getter_optimizer_bool = False
        self._input_sequence_getter_optimizer_last = col_names
        return self._input_sequence
