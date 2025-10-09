from importlib.resources import as_file, files
from os.path import join as pj
from pathlib import Path

from .cluster import *
from .comparison import *
from .compound import *
from .convolution import *
from .data_transform import *
from .distribution import *
from .enums import *
from .error import *
from .estimate import *
from .histogram import *
from .mixture import *
from .output import *
from .plot import *
from .regression import *
from .simulate import *
from .vectors import *


def get_shared_data(file):
    import openalea.stat_tool

    datadir = files("openalea.stat_tool.data")
    with as_file(datadir / file) as f:
        return str(f)
