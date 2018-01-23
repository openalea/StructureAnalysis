import autowig
import os
import shutil
import sys
try:
    from path import Path
except:
    from path import path as Path


prefix = Path(sys.prefix).abspath()
root = Path('.')

os.makedirs(root/'include'/'stat_tool')
for file in (root/'src'/'cpp').walkfiles('*.h*'):
    file.copy(root/'include'/'stat_tool'/file.name)

headers = list((root/'include'/'stat_tool').walkfiles('*.h'))

flags = ['-x', 'c++', '-std=c++11']
flags.append('-I' + str((prefix/'include').abspath()))
flags.append('-I' + str((root/'include').abspath()))

asg = autowig.AbstractSemanticGraph()    
autowig.parser.plugin = 'clanglite'
asg = autowig.parser(asg, headers,
                          flags = flags,
                          bootstrap = 2)

def stat_tool_controller(asg):
    for noncopyable in ['class ::std::basic_streambuf< char, struct ::std::char_traits< char > >',
                        'class ::std::codecvt< char, char, __mbstate_t >',
                        'class ::std::locale::facet',
                        'class ::std::locale::id',
                        'class ::std::ctype< char >',
                        'class ::std::ios_base',
                        'class ::std::basic_istream< char, struct ::std::char_traits< char > >',
                        'class ::std::basic_ostream< char, struct ::std::char_traits< char > >',
                        'class ::std::basic_ostringstream< char, struct ::std::char_traits< char >, class ::std::allocator< char > >',
                        'class ::std::basic_ios< char, struct ::std::char_traits< char > >',
                        'class ::std::basic_stringbuf< char, struct ::std::char_traits< char >, class ::std::allocator< char > >']:
        asg[noncopyable].is_copyable = False
    for cls in asg.classes():
        for fld in cls.fields(access='public'):
            if fld.qualified_type.unqualified_type.globalname == 'class ::std::locale::id':
                fld.boost_python_export = False
    for specialization in asg['class ::std::reverse_iterator'].specializations():
        specialization.boost_python_export = False
    asg['::std::ios_base::openmode'].qualified_type.boost_python_export = True
    return asg

autowig.controller['stat_tool'] = stat_tool_controller
autowig.controller.plugin = 'stat_tool'
asg = autowig.controller(asg)

autowig.generator.plugin = 'boost_python_internal'
wrappers = autowig.generator(asg,
                             module = root/'src'/'py'/'wrapper'/'_stat_tool.cpp',
                             decorator = root/'src'/'py'/'stat_tool'/'_stat_tool.py',
                             prefix = 'wrapper_')

wrappers.write()

shutil.rmtree('include')

import pickle
with open('ASG.pkl', 'w') as filehandler:
    pickle.dump(asg, filehandler)