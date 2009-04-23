"""setup file for tree reduction package"""
import os, sys
from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
from os.path import join as pj
__revision__ = "$Id: setup.py 6086 2009-03-13 16:24:30Z cokelaer $"



name = 'VPlants.tree_reduction'
namespace = 'openalea'
pkg_name = 'openalea.tree_reduction'
version = '0.6.2'
description = 'tree reduction.'
long_description = """"""
author = 'Farah Ben-naoum'
author_email = ''
url = 'http://openalea.gforge.inria.fr'
license = 'Cecill-C'

packages = [ namespace+"."+pkg for pkg in find_packages('src') if 'openalea' not in pkg]

setup(
    name=name,
    version=version,
    description=description,
    author=author,
    author_email=author_email,
    url=url,
    license=license,

    namespace_packages=['openalea'],
    create_namespaces = True,
    zip_safe = False,

    packages = packages,
    package_dir={ pkg_name : pj('src','tree_reduction'), '' : 'src' },

    # Dependencies
    install_requires = ['openalea.deploy'],
    dependency_links = ['http://openalea.gforge.inria.fr/pi'],

    )


    
