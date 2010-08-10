import os.path
from os.path import join as pj

cdir = os.path.dirname(__file__)
pdir = pj(cdir, "..", "..", "stat_tool")
pdir = os.path.abspath(pdir)
import openalea.stat_tool
from openalea.stat_tool import __path__
__path__ = [pdir] + __path__[:]

from openalea.stat_tool.__init__ import *


"""
# Redirect path
import os
import vplants.stat_tool
__path__ = vplants.stat_tool.__path__ + __path__[:]

from vplants.stat_tool.__init__ import *
"""
