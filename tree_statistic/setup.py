# -*- coding: utf-8 -*-
__revision__ = "$Id$"

import os, sys
from os.path import join as pj

from setuptools import setup, find_packages

from openalea.deploy.metainfo import read_metainfo

metadata = read_metainfo('metainfo.ini', verbose=True)
for key,value in metadata.iteritems():
    exec("%s = '%s'" % (key, value))


# Scons build directory
build_prefix = "build-scons"
scons_parameters=["build_prefix="+build_prefix]

# platform dependencies
install_requires = ['vplants.aml','vplants.tree']
setup_requires = install_requires + ['openalea.deploy']
          
packages=[namespace,
          namespace+".tree_statistic",
          namespace+".tree_statistic.trees",
          namespace+".tree_statistic.hmt",
          namespace+".tree_statistic.int_fl_containers",
          #namespace+".treestat_wralea",
         ]

if __name__ == '__main__':
    
    setup(name=name,
          version=version,
          author=authors,
          description=description,
          url=url,
          license=license,
          
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
        
          # Packages
          #namespace_packages = [namespace],
          #create_namespaces = False,

          packages=packages,


          package_dir={
            '' : 'src',
            },

          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so', '*.dylib'],
                           namespace+".tree_statistic": ['*/*.pyd', '*/*.so', '*.dylib']},
          zip_safe=False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = {'include' : pj(build_prefix, 'include') },

          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          
          #entry_points = {
          #  "wralea": ["tree_statistic = openalea.treestat_wralea",
          #              "macro = openalea.treestat_wralea.macro",
          #              "demo = openalea.treestat_wralea.demo",
          #             ]
          #  },
          pylint_packages = ['src'+os.sep+package.replace('.', os.sep) for package in packages]
            )
    
