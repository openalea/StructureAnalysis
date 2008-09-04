# Redirect path
import os
import openalea.stat_tool

__path__ = openalea.stat_tool.__path__ + __path__[:]
from openalea.stat_tool.__init__ import *
