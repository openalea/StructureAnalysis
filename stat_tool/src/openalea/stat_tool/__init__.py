# Redirect path
import os

cdir = os.path.dirname(__file__)
pdir = os.path.join(cdir, "../../stat_tool")
pdir = os.path.abspath(pdir)

__path__ = [pdir] + __path__[:]

from vplants.stat_tool.__init__ import *

