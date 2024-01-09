"""Testing class EyeMovementData."""
from ema import EyeMovementData
import pytest


@pytest.fixture
def eye_movement_data_init():
    return EyeMovementData()


def test_remove_scanpaths_shorter_than_x_fixations(eye_movement_data_init):
    x = 8
    eye_movement_data_init.remove_scanpaths_shorter_than_x_fixations(x)
    gb = eye_movement_data_init.dataframe.groupby(['SUBJ', 'TEXT'], as_index=False).size() < x
    assert sum(gb) == 0
