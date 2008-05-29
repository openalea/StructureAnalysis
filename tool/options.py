
#scons parameters file
#use this file to pass custom parameter to SConstruct script

import sys
if('win' in sys.platform):
    compiler='mingw'

else:
    QTDIR="/usr"
    QT_VERSION=4
