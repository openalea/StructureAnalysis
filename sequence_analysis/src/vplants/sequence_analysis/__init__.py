# Redirect path
import os
import openalea.sequence_analysis

__path__ = openalea.sequence_analysis.__path__ + __path__[:]
from openalea.sequence_analysis.__init__ import *
