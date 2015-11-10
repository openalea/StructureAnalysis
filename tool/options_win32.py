#scons parameters file
#use this file to pass custom parameter to SConstruct script

import sys
if('win' in sys.platform):
    compiler='mingw'
    if compiler == 'mingw':
        gl_includes = r'C:\MinGW\include\GL'
        gl_lib = r'C:\MinGW\lib' 
        QTDIR=r'C:\Qt\4.2.3'
        QT_VERSION='4'
    else: 
        EXTRA_LIBS = 'advapi32.lib user32.lib gdi32 shell32.lib'

