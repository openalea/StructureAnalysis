# -*- coding: utf-8 -*-
__revision__ = "$Id$"

"""setup file for stat_tool package"""

import os, sys
from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
from openalea.deploy.setup import *
from os.path import join as pj

from openalea.deploy.metainfo import read_metainfo
metadata = read_metainfo('metainfo.ini', verbose=True)
for key,value in metadata.iteritems():
    exec("%s = '%s'" % (key, value))


build_prefix = "build-scons"

# Scons build directory
scons_parameters = [] #["build_prefix=" + build_prefix]


# platform dependencies
install_requires = []#[binary_deps('vplants.tool')]
if sys.platform.startswith('win'):
    install_requires += [binary_deps("boost")]
install_requires = []
setup_requires = install_requires + ['openalea.deploy']

OPENALEA_NAMESPACE = False
namespace_packages = ['structure_analysis']
create_namespaces = False
packages = ['structure_analysis',
            'structure_analysis.stat_tool'
           ]

package_dir={}
package_dir[""] = pj("src", "py")
package_dir["structure_analysis"] = pj("src", "py", "structure_analysis")
package_dir["structure_analysis.stat_tool"] = pj("src", "py", "structure_analysis", "stat_tool")

if OPENALEA_NAMESPACE:
    namespace_packages.append('openalea')
    create_namespaces = True

    packages.extend(['openalea', 'openalea.stat_tool'])
    package_dir["openalea.stat_tool"] = pj("src", "stat_tool")
    package_dir[""] = "src"


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

          namespace_packages=namespace_packages,
          create_namespaces=create_namespaces,

          # Packages
          packages=packages,
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
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],


       )



