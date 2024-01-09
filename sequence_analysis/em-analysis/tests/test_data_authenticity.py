# -*- coding: utf-8 -*-

"""
Help testing data authenticity.
File name must be manually passed as an EyeMovementData() argument and handled by the constructor.
It should also respect a certain formatting in terms of column name. (ISLAST, ISFIRST, READMODE, SUBJ, TEXT).
"""

# from ema.config import DATA_PATH
from ema import EyeMovementData
import numpy as np
# import os.path
import pytest


@pytest.fixture(scope="module")
def eye_movement_data_init():
    return EyeMovementData()  # os.path.join(DATA_PATH, 'oculo_data_nov17.xlsx')


def test_is_last(eye_movement_data_init):
    df = eye_movement_data_init.dataframe
    n_samples = df.groupby(['SUBJ', 'TEXT'], as_index=False).size().shape[0]
    n_is_last = (df['ISLAST'] == 1).sum()
    assert n_samples == n_is_last


def test_is_first(eye_movement_data_init):
    df = eye_movement_data_init.dataframe
    n_samples = df.groupby(['SUBJ', 'TEXT'], as_index=False).size().shape[0]
    n_is_first = (df['ISFIRST'] == 1).sum()
    assert n_samples == n_is_first


def test_to_be_replaced(eye_movement_data_init):
    df = eye_movement_data_init.dataframe
    n_samples = df.groupby(['SUBJ', 'TEXT'], as_index=False).size().shape[0]
    n_to_be_replaced = (df['READMODE'].astype(str).str.contains('_TOBEREPLACED_')).sum()
    assert n_to_be_replaced == n_samples or n_to_be_replaced == 0


def test_constant_column(eye_movement_data_init):
    df = eye_movement_data_init.dataframe
    for col in df.columns:
        assert df[col].unique().shape[0] != 1


def test_replicates(eye_movement_data_init):
    df = eye_movement_data_init.dataframe
    assert np.all(np.array(df.groupby(['SUBJ', 'TEXT'])['ISLAST'].sum() == 1))
