# -*- coding: utf-8 -*-

import os, sys
from os.path import join as pj

from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps

name = 'VPlants.Sequence_analysis'
namespace = 'openalea'
version = "0.7.0"
description = 'sequence analysis library',
long_description = """Python version of sequence analysis (AML stat)"""
authors = 'Y. Guédon, T. Cokelaer',
author_email = 'yann.guedon@cirad.fr, christophe.pradal@cirad.fr, thomas.cokelaer@sophia.inria.fr'
url = 'http://www-sop-inria.fr/virtualplants'
license = 'GPL'
__revision__ = "$Id: setup.py 6086 2009-03-13 16:24:30Z cokelaer $"




packagename = 'sequence_analysis'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
install_requires = [binary_deps('vplants.stat_tool')]
#install_requires = []
setup_requires = install_requires + ['openalea.deploy']

if sys.platform.startswith('win'):
    install_requires += [binary_deps("boostpython")]
    setup_requires += [binary_deps("boostpython")]


if __name__ == '__main__':

    setup(name=name,
          version=version,
          author=authors,
          author_email=author_email,
          description=description,
          long_description=long_description,
          url=url,
          license=license,


          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters
          scons_parameters=scons_parameters,

          namespace_packages=['openalea'],
          create_namespaces=True,

          # Packages
          packages=[ 
                    'openalea.sequence_analysis',
                    'vplants'],
          

           #package_dir={'openalea.sequence_analysis' : 'src/openalea/sequence_analysis',
                       #'' : 'build/lib', # hack to use develop command
#                       },

#          package_dir={ "" : "src", vplants':'src/vplants/sequence_analysis' },
          package_dir={ "" : "src",  },


          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so', '*.dylib'],},
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },


          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )

