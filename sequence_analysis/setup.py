# -*- coding: utf-8 -*-
__revision__ = "$Id$"

import os, sys
from os.path import join as pj

from setuptools import setup, find_namespace_packages
#from openalea.deploy.binary_deps import binary_deps


# from openalea.deploy.metainfo import read_metainfo
# metadata = read_metainfo('metainfo.ini', verbose=True)
# for key,value in metadata.iteritems():
#     exec("%s = '%s'" % (key, value))

# Meta-information
name='OpenAlea.SequenceAnalysis'
version='2.0.0'
description='sequence analysis library'
long_description='Python version of sequence analysis (AML stat).'
authors='Y. Guedon, JB. Durand, P. Fernique, C. Pradal, T. Cokelaer'
authors_email='christophe.pradal@cirad.fr'
url='https://github.com/openalea/StructureAnalysis/'
license='CeCILL-C'


build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
install_requires = []
setup_requires = install_requires + ['openalea.deploy']

namespace = 'openalea'
packages = find_namespace_packages(where='src', include=['openalea.*'])
package_dir = {'': 'src'}

"""
if sys.platform.startswith('win'):
    install_requires = [binary_deps('openalea.stattool')]
    install_requires += [binary_deps("boost")]
    setup_requires += [binary_deps("boost")]
"""

if __name__ == '__main__':

    setup(name=name,
          version=version,
          author=authors,
          author_email=authors_email,
          description=description,
          long_description=long_description,
          url=url,
          license=license,


          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters
          scons_parameters=scons_parameters,

          namespace_packages=['openalea'],
          #create_namespaces=True,

          # Packages
          packages=packages,

          package_dir=package_dir,
          share_dirs = { 'share' : 'share' },


          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so', '*.dylib'],},
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },

          entry_points = {
            "wralea": ["openalea.stats = sequence_analysis_wralea",
                       "openalea.demo = sequence_analysis_wralea.demo.change_point",
            ]
            },

          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          #dependency_links = ['http://openalea.gforge.inria.fr/pi'],

          #pylint_packages = ['src/openalea/sequence_analysis']
          )

