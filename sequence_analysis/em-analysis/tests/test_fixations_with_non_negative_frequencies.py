# -*- coding: utf-8 -*-
"""Compare extraction of fixations with negative word frequencies:
Excel way vs pandas"""
__version__ = ""

from ema import MODELS_PATH
from ema import DATA_PATH
from ema import EyeMovementData
from ema import Model
from nose import with_setup
from openalea.sequence_analysis import HiddenSemiMarkov, Sequences, Shift, SelectVariable
from openalea.sequence_analysis._sequence_analysis import _MarkovianSequences
import os


def init():
    return setup_func()


def setup_func():
    global s3a, s3b, model, model_init
    # import data (pandas)
    data = EyeMovementData()
    # remove fixations that have negative word frequencies
    data.drop_negative_wfreq_fixations()
    # data filtered with excel
    excel_data = os.path.join(DATA_PATH, "data.txt")
    f = open(excel_data)
    c = f.read()
    f.close()
    l = eval(c)
    s = Sequences(l)
    s = Shift(s, 3, -1)  
    s = _MarkovianSequences(s)
    s3a = SelectVariable(s, 3)
    # Prepare initialization and extract data from model
    hsmc_file = os.path.join(MODELS_PATH, "model3_init5abs.hsmc")
    model = Model(data, init_hsmc_file=hsmc_file)
    model_init = HiddenSemiMarkov(hsmc_file)
    s3b = model._sequential_data.get_input_sequence
    return s3a, s3b, model, model_init


@with_setup(setup_func)      
def test_nb_sequences():
    """Test whether both data sets have same number of sequences"""   
    msg = "Data sets have different numbers of sequences: "
    msg += "Excel: " + str(len(s3a)) + "; "
    msg += "Pandas: " + str(len(s3b))
    assert len(s3a) == len(s3b), msg


@with_setup(setup_func)
def test_sequence_lengths():
    """Test whether all sequences in both data sets have same lengths"""   
    msg = "At least two sequences have different lengths: ("
    msge = "Excel: "
    msgp = "Pandas: "
    i = 0
    while len(s3a[i]) == len(s3b[i]) and i < min(len(s3a)-1, len(s3b)-1):
        i += 1
    msg += "sequence " + str(i) + ") " 
    msg += msge + str(len(s3a[i])) + "; "
    msg += msgp + str(len(s3b[i]))
    assert i == min(len(s3a)-1, len(s3b)-1), msg


@with_setup(setup_func)
def test_sequence_contents():
    """Test whether all sequences in both data sets have same contents"""   
    msg = "At least two sequences have different contents: ("
    msge = "Excel: "
    msgp = "Pandas: "
    i = 0
    while str(s3a[i]) == str(s3b[i]) and i < min(len(s3a)-1, len(s3b)-1):
        i += 1
    msg += "sequence " + str(i) + ") " 
    msg += msge + str(s3a[i]) + "; "
    msg += msgp + str(s3b[i])
    assert i == min(len(s3a)-1, len(s3b)-1), msg


if __name__ == "__main__":
    s3a, s3b, model, model_init = init()
    test_nb_sequences()
    test_sequence_lengths()
    test_sequence_contents()

