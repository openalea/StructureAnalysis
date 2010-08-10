# Redirect path
import os

cdir = os.path.dirname(__file__)
pdir = os.path.join(cdir, "../../sequence_analysis")
pdir = os.path.abspath(pdir)

__path__ = [pdir] + __path__[:]

from vplants.sequence_analysis.__init__ import *

