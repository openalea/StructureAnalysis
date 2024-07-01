# -*- coding: utf-8 -*-
"""setup file for stat_tool package"""

import os, sys
from setuptools import setup, find_namespace_packages
#from openalea.deploy.binary_deps import binary_deps
from openalea.deploy.setup import *
from os.path import join as pj

# from openalea.deploy.metainfo import read_metainfo
# metadata = read_metainfo('metainfo.ini', verbose=True)
# for key,value in metadata.iteritems():
#     exec("%s = '%s'" % (key, value))

# Meta-information
name='OpenAlea.StatTool'
version='2.0.0'
description='Basic Statistical tools'
long_description='Basic Statistical tools used by different Structure Analysis libraries.'
authors='Y. Guedon, JB. Durand, P. Fernique, C. Pradal, T. Cokelaer, S. Dufour'
authors_email='christophe.pradal@cirad.fr'
url='https://github.com/openalea/StructureAnalysis/'
license='CeCILL-C'

build_prefix = "build-scons"

# Scons build directory
scons_parameters = ["build_prefix=" + build_prefix]


# platform dependencies
install_requires = []
setup_requires = ['openalea.deploy']

namespace = 'openalea'
packages = find_namespace_packages(where='src', include=['openalea.*'])
package_dir = {'': 'src'}


if __name__ == '__main__':

    setup(name=name,
          version=version,
          description=description,
          long_description=long_description,
          author=authors,
          author_email=authors_email,
          url=url,
          license=license,
          platforms = platforms,

          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  v
          scons_parameters=scons_parameters,

          namespace_packages=[namespace], #, "structure_analysis"],
          #namespace_packages=["structure_analysis"],
          #create_namespaces=False,

          # Packages
          packages=packages,
                    #'structure_analysis',
                    #'structure_analysis.stat_tool',
                    #'stat_tool'
                    #],

          package_dir=package_dir,
          share_dirs = { 'share' : 'share' },


          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so', '*.dylib', '*.png', '*.hsc', '*.seq', '*.aml'],},

          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },


          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          #dependency_links = ['http://openalea.gforge.inria.fr/pi'],


       )



