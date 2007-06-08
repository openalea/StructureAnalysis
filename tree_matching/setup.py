from openalea.distx import setup, find_packages, find_package_dir, Shortcut 
from os.path import join as pj
 
packagename = 'tree_matching'
namespace = 'openalea'

# Scons build directory
build_prefix= "build-scons"

if __name__ == '__main__':
    
    setup(name=packagename,
          version='0.1',
          author='',
          author_email='',
          description='tree matching library',
          url='',
          license='GPL',
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=["build_prefix="+build_prefix],
        
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

