from vplants.autowig import autowig

from path import path

rootdir = path('.').abspath()

includedir = rootdir + "/stat_tool/build-scons/include"
headers = [str(f) for f in includedir.walkfiles('*.h')]

flags = ['-x', 'c++', '-g', '-std=c++11', '-stdlib=libstdc++',
         '-I/usr/include', '-I' + str(rootdir + '/tool/build-scons/include'), '-I' + str(includedir), '-D__STDC_CONSTANT_MACROS',
         '-I/usr/local/lib/clang/3.7.0/include', '-D__STDC_LIMIT_MACROS']

asg = autowig.AbstractSemanticGraph()
autowig.front_end.plugin = 'pyclanglite'
autowig.front_end(asg, headers, flags=flags, silent=True,
                  force_overload=True, bootstrap=3, cache=rootdir + '/.AutoWIG', force=True)

autowig.middle_end.plugin = 'default'
autowig.middle_end(asg)

#autowig.held_type.plugin = 'ptr'

autowig.back_end.plugin = 'boost_python:std_filter'
autowig.back_end(asg)

autowig.back_end.plugin = 'boost_python:export'
autowig.back_end(asg,
                 directory = rootdir + '/stat_tool/src/wrapper',
                 pattern = '.*stat_tool.*',
                 prefix = '_')

#autowig.back_end.plugin = 'boost_python:closure'
#autowig.back_end(asg)
#
#autowig.back_end.plugin = 'boost_python:export'
#autowig.back_end(asg,
#                 directory = rootdir + '/misc/src/wrapper',
#                 prefix = '_')

autowig.back_end.plugin = 'boost_python:module'
autowig.back_end(asg,
                 filename = rootdir + '/stat_tool/src/wrapper/__stat_tool.cpp',
                 package = 'openalea.stat_tool')
#autowig.back_end(asg,
#                 filename = rootdir + '/misc/src/wrapper/__misc.cpp',
#                 package = 'structure_analysis')

autowig.back_end.plugin = 'boost_python:import'
autowig.back_end(asg,
                 filename = rootdir + '/stat_tool/src/openalea/stat_tool/_stat_tool.py',
                 module = rootdir + '/stat_tool/src/wrapper/__stat_tool.cpp')
#autowig.back_end(asg,
#                 filename = rootdir + '/misc/src/structure_analysos/_misc.py',
#                 module = rootdir + '/misc/src/wrapper/__misc.cpp')

#autowig.back_end.plugin = 'cpp:char_ptr_to_string'
#autowig.back_end(asg, '^::stat_tool::.*', rootdir + '/stat_tool/src/cpp/char_ptr_to_string')

autowig.back_end.plugin = 'on_disk'
autowig.back_end(asg,
                 pattern = rootdir + '/(stat_tool|misc)/src/(wrapper|openalea)/.*')
