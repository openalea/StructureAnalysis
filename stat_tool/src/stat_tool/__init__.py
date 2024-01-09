from distribution import *
from histogram import *
from mixture import *
from multivariate_mixture import *
from compound import *
from convolution import *
from vectors import *
from estimate import *
from simulate import *
from comparison import *
from regression import *
from output import *
from data_transform import *
from cluster import *
from enums import *
from plot import *
from error import *

from openalea.deploy.shared_data import get_shared_data_path
from os.path import join as pj

def get_shared_data(file):
    import openalea.stat_tool
    shared_data_path = get_shared_data_path(openalea.stat_tool.__path__)
    return pj(shared_data_path, file)
