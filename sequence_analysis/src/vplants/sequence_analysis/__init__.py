import os.path
from os.path import join as pj

cdir = os.path.dirname(__file__)
pdir = pj(cdir, "..", "..", "sequence_analysis")
pdir = os.path.abspath(pdir)
import openalea.sequence_analysis
from openalea.sequence_analysis import __path__
__path__ = [pdir] + __path__[:]

from openalea.sequence_analysis.__init__ import *


