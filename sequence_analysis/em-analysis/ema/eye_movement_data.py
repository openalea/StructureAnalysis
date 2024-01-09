"""Stores oculometric data as pandas dataframe."""
from .config import (IS_FIRST_COL, IS_LAST_COL, OFF_DURATION_COL, SUBJECT_NAME_COL, TEXT_NAME_COL,
                     CHAR_INCREMENT_COL, WORD_INCREMENT_COL, SACCADE_DIRECTION_COL,
                     SACCADE_AMPLITUDE_COL, FIXATION_DURATION_COL, READMODE_COL, WORD_FREQUENCY_COL)
import math
from openalea.sequence_analysis import Sequences
from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences
import numpy as np
import pandas as pd

"""
TODO:
* Texts A cost inst * wfreq
* cmp pd.DataFrame, noseq, nofix, segmentation (hist phases)
"""

class EyeMovementDataFrame(pd.DataFrame):
    def __init__(self, *args, **kwargs):
        index_as_ms = kwargs.pop('index_as_ms', False)
        super(EyeMovementDataFrame, self).__init__(*args, **kwargs)
        self.index_as_ms = index_as_ms if type(index_as_ms) == bool else None
        self._input_sequence = None
        self._input_sequence_getter_optimizer_last = None
        self._input_sequence_getter_optimizer_bool = False
        self.restored_data = None

    def col_to_seq(self, col_names):
        """Extract cols from DataFrame as a sequence of sequence."""
        if all(isinstance(item, str) for item in col_names):
            for col_name in col_names:
                assert(col_name in self.columns)
        else:
            raise ValueError('All %s must be in %s' % (col_names, self.column_names))
        data = []
        is_last = self[IS_LAST_COL] == 1
        scanpath_start_index = 0
        for i in self.index[is_last]:
            data.append(self.loc[range(scanpath_start_index, i+1), col_names].values.tolist())
            scanpath_start_index = i+1
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

    def test_data_authentification(self):
        print('ISFIRST: ', (self[IS_FIRST_COL] == 1).sum())
        print('ISLAST: ', (self[IS_LAST_COL] == 1).sum())
        print('SUBJ_NAME/TEXT: ', self.groupby([SUBJECT_NAME_COL, TEXT_NAME_COL]).count().shape[0])
        assert (self[IS_FIRST_COL] == 1).sum() == (self[IS_LAST_COL] == 1).sum() == \
               self.groupby([SUBJECT_NAME_COL, TEXT_NAME_COL]).count().shape[0]

    def drop_pathological_fixations(self):
        df = self.copy()
        set_is_last = []
        # in case the fixations before the last one also have wfreq < 0
        for i in df.index[(df.OFF_DUR == 1) & (df.ISLAST == 1)]:
            j = i - 1
            while df.at[j, 'OFF_DUR'] == 1:
                j -= 1
            set_is_last.append(j)

        set_is_first = []
        for i in df.index[(df.OFF_DUR == 1) & (df.ISFIRST == 1)]:
            j = i + 1
            while df.at[j, 'OFF_DUR'] == 1:
                j += 1
            set_is_first.append(j)

        df.at[set_is_last, 'ISLAST'] = 1
        df.at[set_is_first, 'ISFIRST'] = 1
        df = df[df.OFF_DUR == 0]
        df = df.reset_index(drop=True)
        return df

    def remove_scanpaths_with_more_than_x_missing_values(self, x):
        df = self.copy()
        gb = df.groupby(['SUBJ_NAME', 'TEXT'])['CINC'].apply(lambda x: pd.isnull(x).sum()) >= x
        for i in gb.index:
            if gb[i]:
                idx_to_drop = df[(df['SUBJ_NAME'] == i[0]) & (df['TEXT'] == i[1])].index
                df = df.drop(idx_to_drop)
        df = df.reset_index(drop=True)
        return df

    def compute_saccade_direction(self):
        """Compute saccade direction."""
        self['SACDIR'] = 'forward'
        for i in self.index:
            if self.at[i, 'ISLAST'] == 0:
                diff_x = self.at[i + 1, 'X'] - self.at[i, 'X']
                diff_y = self.at[i + 1, 'Y'] - self.at[i, 'Y']
                if diff_y > 30:
                    self.at[i, 'SACDIR'] = 'downward'
                elif diff_y < -30:
                    self.at[i, 'SACDIR'] = 'upward'
                elif diff_x < 0:
                    self.at[i, 'SACDIR'] = 'backward'
            else:
                self.at[i, 'SACDIR'] = 'last'

    def recompute_saccade_amplitude(self):
        """Set saccade amplitude to np.nan for line breaks and last outgoing saccade"""
        self['SACAMP'] = abs(self['SACAMP'])
        for i in self.index:
            if self.at[i, 'ISLAST'] == 1:
                self.at[i, 'SACAMP'] = np.nan
            elif self.at[i + 1, 'ISLAST'] != 1:
                i += 1
                diff_y = self.at[i + 1, 'Y'] - self.at[i, 'Y']
                a = math.sqrt((self.at[i + 1, 'X'] - self.at[i, 'X']) ** 2
                              + (self.at[i + 1, 'Y'] - self.at[i, 'Y']) ** 2)
                b = math.sqrt((self.at[i + 1, 'X'] - self.at[i - 1, 'X']) ** 2
                              + (self.at[i + 1, 'Y'] - self.at[i - 1, 'Y']) ** 2)
                c = math.sqrt((self.at[i - 1, 'X'] - self.at[i, 'X']) ** 2
                              + (self.at[i - 1, 'Y'] - self.at[i, 'Y']) ** 2)
                angle = math.acos((a ** 2 + c ** 2 - b ** 2) / (2 * a * c)) * (180. / math.pi)
                if diff_y > 30 and angle < 20:  # go to next line
                    self.at[i, 'SACAMP'] = np.nan

    def drop_negative_wfreq_fixations(self):
        df = self.copy()
        is_last = df['ISLAST'] == 1
        neg_wfreq = df['WFREQ'] < 0
        set_is_last = []
        # in case the fixations before the last one also have wfreq < 0
        for i in df.index[is_last & neg_wfreq]:
            j = i-1
            while df.loc[j, 'WFREQ'] < 0:
                j -= 1
            set_is_last.append(j)

        df.loc[set_is_last, 'ISLAST'] = 1
        df = self.drop(self.index[neg_wfreq])
        df = self.reset_index(drop=True)
        df = self.drop(9947)  # wtf ?
        df = self.reset_index(drop=True)
        return df

    def remove_scanpaths_shorter_than_x_fixations(self, x):
        df = self.copy()
        gb = df.groupby(['SUBJ_NAME', 'TEXT'], as_index=False).size() < x
        for i in gb.index:
            if gb[i]:
                idx_to_drop = df[(df['SUBJ_NAME'] == i[0]) & (df['TEXT'] == i[1])].index
                df = df.drop(idx_to_drop)
        df = df.reset_index(drop=True)
        return df
