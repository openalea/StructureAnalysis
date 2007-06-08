from openalea.distx import setup, find_packages, find_package_dir, Shortcut 
from os.path import join as pj
import os
 
packagename = 'tree_statistic'
namespace = 'openalea'

# Scons build directory
build_prefix= "build-scons"
scons_parameters=["build_prefix="+build_prefix]
if os.name == 'nt':
    scons_parameters=["build_prefix="+build_prefix,"compiler=mingw"]
    

if __name__ == '__main__':
    
    setup(name=packagename,
          version='0.1',
          author='J.-B. Durand',
          author_email='jean-baptiste.durand@cirad.fr',
          description='tree statistic',
          url='',
          license='GPL',
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
 
        
          # Packages
          packages=find_packages(where='src', namespace=namespace),
          package_dir=find_package_dir(where='src', namespace=namespace), 
      
          # Add package platform libraries if any
          include_package_lib=True,
                    
	  # Copy shared data in default OpenAlea directory
	  # Map of 'destination subdirectory' : 'source subdirectory'
	  external_data={pj('include'):  pj(build_prefix, 'include'),
                         pj('lib'):  pj(build_prefix, 'lib'),
                         },
 
)
