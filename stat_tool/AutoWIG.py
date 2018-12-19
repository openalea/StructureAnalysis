import autowig
import os
import shutil
import sys

from path import Path


prefix = Path(sys.prefix).abspath()
root = Path('.')

inc_dir = root/'include'/'stat_tool'
if not inc_dir.exists():
    os.makedirs(root/'include'/'stat_tool')
else:
    for file in inc_dir.files('*.h*'):
        file.remove()

for file in (root/'src'/'cpp').walkfiles('*.h*'):
    file.copy(root/'include'/'stat_tool'/file.name)

headers = list(inc_dir.walkfiles('*.h'))

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
        if noncopyable in asg:
            asg[noncopyable].is_copyable = False
    for cls in asg.classes():
        for fld in cls.fields(access='public'):
            if fld.qualified_type.unqualified_type.globalname == 'class ::std::locale::id':
                fld.boost_python_export = False
    for specialization in asg['class ::std::reverse_iterator'].specializations():
        specialization.boost_python_export = False
    for ctr in asg['class ::std::basic_string< char, struct ::std::char_traits< char >, class ::std::allocator< char > >'].constructors():
        ctr.boost_python_export = False
    asg['::std::ios_base::openmode'].qualified_type.boost_python_export = True


    for cls in asg['class ::std::vector'].specializations(partial = False):
        for method in cls.methods():
            if method.localname in ['resize', 'shrink_to_fit', 'operator[]']:
                if isinstance(method.boost_python_export, bool):
                    method.boost_python_export = False
        for constructor in cls.constructors():
            if not(constructor.nb_parameters == 0 or constructor.nb_parameters == 1 and constructor.parameters[0].qualified_type.unqualified_type == cls):
                if isinstance(constructor.boost_python_export, bool):
                    constructor.boost_python_export = False

    for cls in asg['class ::std::allocator'].specializations(partial = False):
        cls.boost_python_export = False
    if 'class ::std::reverse_iterator' in asg:
        for cls in asg['class ::std::reverse_iterator'].specializations(partial = False):
            cls.boost_python_export = False
    if 'class ::std::initializer_list' in asg:
        for cls in asg['class ::std::initializer_list'].specializations(partial = False):
            cls.boost_python_export = False
    if 'class ::std::default_delete' in asg:
        for cls in asg['class ::std::default_delete'].specializations(partial = False):
            cls.boost_python_export = False
    for mtd in asg['::std::string'].qualified_type.desugared_type.unqualified_type.methods():
        if mtd.localname in ['substr', 'compare']:
            mtd.boost_python_export = False


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

# strange line
#shutil.rmtree('include')

import pickle
with open('ASG.pkl', 'w') as filehandler:
    pickle.dump(asg, filehandler)
