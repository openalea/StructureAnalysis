"""setup file for stat_tool package"""
import os, sys
from setuptools import setup, find_packages
from openalea.deploy.binary_deps import binary_deps
from os.path import join as pj



name = 'VPlants.Stat_Tool'
namespace = 'openalea'
# to get the version
version="0.7.0"
description = 'statistics' 
long_description = """todo"""
author = 'Y. Guedon, JB. Durand, T. Cokelaer',
author_email = 'samuel.dufour@sophia.inria.fr, christophe.pradal@cirad.fr, thomas.cokelaer@sophia.inria.fr'
url = 'http://openalea.gforge.inria.fr/doc/vplants/stat_tool/doc/html/contents.html'
license = 'GPL' 
__revision__ = "$Id$"


packagename = 'stat_tool'
build_prefix = "build-scons"

# Scons build directory
scons_parameters = ["build_prefix=" + build_prefix]


# platform dependencies
install_requires = [binary_deps('vplants.tool')]
if sys.platform.startswith('win'):
    install_requires += [binary_deps("boostpython")]

setup_requires = install_requires + ['openalea.deploy']


if __name__ == '__main__':
    
    setup(name=name,
          version=version,
          description=description,
          long_description=long_description,
          author=author,
          author_email=author_email,
          url=url,
          license=license,
 
          
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  v
          scons_parameters=scons_parameters,

          namespace_packages=['openalea'],
          create_namespaces=True,

          # Packages
          packages=['openalea.stat_tool', 
                    'stat_tool_wralea.demo.change_point', 
                    'stat_tool_wralea.demo.stat_tool_tutorial', 
                    'stat_tool_wralea', 
                    'vplants',
                    ],

          package_dir={ "" : "src"  },
          
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
   
           # entry_points
           entry_points = {
            "wralea": ["vplants.stats = stat_tool_wralea",
                       "vplants.demo = stat_tool_wralea.demo.change_point",
            ]
            },

#            pylint_packages = ['src/openalea/stat_tool'],
#            pylint_options = ['--disable-msg=C0103']

       )


    
