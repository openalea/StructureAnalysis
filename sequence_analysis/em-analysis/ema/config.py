# -*- coding: utf-8 -*-
__author__ = 'Brice Olivier'


import os
cdir = os.path.dirname(__file__)

# OUTPUT_PATH = '/home/bolivier/PyENE/PyENE/em-analysis/output/'
DATA_PATH = os.path.join(cdir, 'share','data')
OUTPUT_PATH = os.path.join(cdir, 'share','models')
MODEL_PATH = OUTPUT_PATH

IS_FIRST_COL = 'ISFIRST'
IS_LAST_COL = 'ISLAST'
OFF_DURATION_COL = 'OFF_DUR'
SUBJECT_NAME_COL = 'SUBJ_NAME'
TEXT_NAME_COL = 'TEXT'
CHAR_INCREMENT_COL = 'CINC'
WORD_INCREMENT_COL = 'WINC'
SACCADE_DIRECTION_COL = 'SACDIR'
SACCADE_AMPLITUDE_COL = 'SACAMP'
FIXATION_DURATION_COL = 'FDUR'
READMODE_COL = 'READMODE'
WORD_FREQUENCY_COL = 'WFREQ'
