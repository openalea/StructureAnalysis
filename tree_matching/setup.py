# -*- coding: utf-8 -*-
__revision__ = "$Id$"

from setuptools import setup
#import os, sys
#from os.path import join as pj
 
packagename = 'tree_matching'
namespace = 'openalea'

    
setup(name="OpenAlea.TreeMatching",
      version='1.0.0',
      author='Pascal Ferraro, Aida Ouangraoua',
      description='Tree matching library',
      url='https://github.com/openalea/StructureAnalysis',
      license='GPL',
                  
          # Packages
          #namespace_packages = ["openalea"],
          #create_namespaces = True,
          
      # pure python  packages
      packages= [ "openalea", 
                  "openalea.tree_matching"
                      ],
          
      # python packages directory
      package_dir= { '' : 'src',
                         },

      # Add package platform libraries if any
      include_package_data=True,
      package_data = {'' : ['*.pyd', '*.so', '*.dylib'],},
      zip_safe=False,

      # Specific options of openalea.deploy
      #lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
      #inc_dirs = {'include' : pj(build_prefix, 'include') },
          
)


    
