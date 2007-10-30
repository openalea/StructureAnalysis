from setuptools import setup, find_packages
import os, sys
from os.path import join as pj
 
packagename = 'tree_statistic'
namespace = 'openalea'
build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]

# platform dependencies
install_requires = ['vplants.aml','vplants.tree']
setup_requires = install_requires + ['openalea.deploy']

if __name__ == '__main__':
    
    setup(name='VPlants.Tree_Statistic',
          version='0.1.1',
          author='JB durand',
          description='Tree statistic library',
          url='',
          license='GPL',
          
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
        
          # Packages
          namespace_packages = [namespace],
          create_namespaces = True,

          packages=[namespace+".stat_tool", namespace+".tree_statistic",],
          package_dir={namespace+".stat_tool" : 'src/stat_tool',
                       namespace+".tree_statistic": 'src/tree_statistic'},
      
          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so'],
                           namespace+".tree_statistic": ['*/*.pyd', '*/*.so']},
          zip_safe=False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = {'include' : pj(build_prefix, 'include') },

          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )


    
