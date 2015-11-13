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
                  force_overload=True, bootstrap=3, cache=rootdir + '/.AutoWIG')

autowig.middle_end.plugin = 'default'
autowig.middle_end(asg)

autowig.back_end.plugin = 'cpp:char_ptr_to_string'
autowig.back_end(asg, '^::stat_tool::.*', rootdir + '/stat_tool/src/cpp/char_ptr_to_string',
        namespace='stat_tool')

autowig.back_end.plugin = 'on_disk'
autowig.back_end(asg,
                 pattern = rootdir + '/stat_tool/src/cpp/char_ptr_to_string.*')
