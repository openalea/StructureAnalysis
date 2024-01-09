# -*- coding: utf-8 -*-
__revision__ = "$Id: setup.py 9869 2010-11-04 17:31:41Z dbarbeau $"

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


# dependencies
setup_requires   = ['openalea.deploy']
install_requires = setup_requires[:] #copy the content not the reference.
if sys.platform.startswith('win'):
    setup_requires += [] 
    #install_requires += ['qt4']

if __name__ == '__main__':
    
    setup(name=name,
          version=version,
          author=authors,
          author_email=authors_email,
          description=description,
          url=url,
          license=license,
 
          # Define where to execute scons
          scons_scripts=['SConstruct'],
          # Scons parameters  
          scons_parameters=scons_parameters,
        
          # Packages
          # packages=
          # package_dir=
      
          # Add package platform libraries if any
          include_package_data=True,
          zip_safe = False,

          # Specific options of openalea.deploy
          lib_dirs = {'lib' : pj(build_prefix, 'lib'),},
          inc_dirs = { 'include' : pj(build_prefix, 'include') },
          

          # Dependencies
          setup_requires   = setup_requires ,
          install_requires = install_requires,
          dependency_links = ['http://openalea.gforge.inria.fr/pi'],
          )

