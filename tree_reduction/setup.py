"""setup file for stat_tool package"""
import os, sys
from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
from os.path import join as pj



name = 'VPlants.tree_reduction'
namespace = 'vplants'
# to get the version
version="0.6.2"
description = 'tree reduction' 
long_description = """"""
author = 'Farah Ben-naoum',
author_email = ''
url = 'http://openalea.gforge.inria.fr'
license = 'GPL' 
__revision__ = "$Id: setup.py 6086 2009-03-13 16:24:30Z cokelaer $"


packagename = 'tree_reduction'

# platform dependencies

setup_requires = ['openalea.core']


if __name__ == '__main__':
    
    setup(name=name,
          version=version,
          description=description,
          long_description=long_description,
          author=author,
          author_email=author_email,
          url=url,
          license=license,
 
#          namespace_packages=['vplants'],
          create_namespaces=False,

          # Packages
          packages=find_packages('src'),

          package_dir={ "" : "src"  },
          
          # Add package platform libraries if any
          include_package_data=True,
          
          zip_safe = False,

          # Dependencies
          setup_requires = setup_requires,
          install_requires = [],
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )


    
