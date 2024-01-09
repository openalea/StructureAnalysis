# This file allow to use openalea packages without installing them.

import os
from os.path import join as pj
from os.path import exists as exists
from distutils.sysconfig import get_python_lib

def set_mod_path(namedir):
    ok = False
    for p in __path__:
        mod_dir = pj(p,namedir)
        if exists(mod_dir):
            ok = True
            break
    if not ok:
        raise "Cannot find "+namedir+" path !"
    bindir = pj(mod_dir,'bin')
    if exists(bindir):
        if not bindir in os.environ['PATH']:
            os.environ['PATH'] = bindir +';' + os.environ['PATH']
    else:
        bindir = pj(mod_dir,'build-scons','bin')
        if exists(bindir):
            if not bindir in os.environ['PATH']:
                os.environ['PATH'] = bindir +';' + os.environ['PATH']
    
    libdir = pj(mod_dir,'lib')
    if os.name == 'posix':
        libpathname = 'LD_LIBRARY_PATH'
    else:
        libpathname = 'PATH'
    if exists(libdir):
        if not libdir in os.environ[libpathname]:
            os.environ[libpathname] = libdir +';' + os.environ[libpathname]
    else:
        libdir = pj(mod_dir,'build-scons','lib')
        if exists(libdir):
            if not os.environ.has_key(libpathname):
                os.environ[libpathname] = libdir
            if not libdir in os.environ[libpathname]:
                os.environ[libpathname] = libdir +';' + os.environ[libpathname]
    return mod_dir
    

def set_pgl_path(namedir = 'PlantGL'):
    pgl_dir = set_mod_path(namedir)
    os.environ['PLANTGLDIR'] = pgl_dir

def set_lsys_path(namedir = 'lpy'):
    set_mod_path(namedir)
    
#set_pgl_path()
#set_lsys_path()
#set_mod_path('fractalysis')

pkg_dirs = [ "PlantGL/src",
             "lpy/src/openalea",
             "celltissue/src",
             "container/grid/src",
             "container/topomesh/src",
             "fractalysis/src",
             "pglviewer/src",
             "svgdraw/src",
             "physics/src"
             ]
#"container/graph/src",

 
def set_path(mod_path):
    c_path = list(mod_path)
    for subdir in  pkg_dirs:
        for p in c_path:
            newpath = os.path.abspath(pj(p, os.path.normpath(subdir)))
            if(os.path.isdir(newpath)):
                mod_path.append( newpath )
                break
    return mod_path

#__path__ = set_path(__path__)

del set_path
del set_mod_path
del set_pgl_path
del set_lsys_path
del get_python_lib
del exists
del pj
