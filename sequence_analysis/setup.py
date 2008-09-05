from setuptools import setup, find_packages
import os, sys
from os.path import join as pj

packagename = 'sequence_analysis'
namespace = 'openalea'

build_prefix = "build-scons"

# Scons build directory
scons_parameters=["build_prefix="+build_prefix]


# platform dependencies
install_requires = ['vplants.stat_tool']
install_requires = []
setup_requires = install_requires + ['openalea.deploy']

if("win" in sys.platform):
    install_requires += ["boostpython"]
    setup_requires += ["boostpython"]


if __name__ == '__main__':

    setup(name='VPlants.Sequence_analysis',
          version='0.3.1',
          author='Yann Guedon',
          description='sequence analysis library',
          url='http://www-sop.inria.fr/virtualplants/',
          license='GPL',

          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters
          scons_parameters=scons_parameters,

          namespace_packages=['openalea'],
          create_namespaces=True,

          # Packages
          packages=['openalea.sequence_analysis'],
          package_dir={'openalea.sequence_analysis' : 'src/sequence_analysis',
                       #'' : 'build/lib', # hack to use develop command
                       },

          # Add package platform libraries if any
          include_package_data=True,
          package_data = {'' : ['*.pyd', '*.so'],},
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },


          # Dependencies
          setup_requires = setup_requires,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )

