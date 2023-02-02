# rm src/py/wrapper/*.cpp
# rm src/py/wrapper/*.h
# source activate vplants && scons; source deactivate vplants
# cd src/cpp && source activate statiskit-dev && python AutoWIG.py; source deactivate statiskit-dev cd ../..
# scons

import autowig
from path import Path
import sys

def controller(asg):
    autowig.controller.plugin = 'default'
    asg = autowig.controller(asg)
    # for function in asg['::statiskit::stl'].functions():
    #     if function.localname in ['generator', 'insert']:
    #         parameter = function.parameters[0].qualified_type.desugared_type
    #         if parameter.is_class:
    #             function.parent = parameter.unqualified_type
    for cls in asg['class ::std::vector'].specializations(partial = False):
        for method in cls.methods():
            if method.localname in ['resize', 'shrink_to_fit', 'operator[]']:
                if isinstance(method.boost_python_export, bool):
                    method.boost_python_export = False
        for constructor in cls.constructors():
            if not(constructor.nb_parameters == 0 or constructor.nb_parameters == 1 and constructor.parameters[0].qualified_type.unqualified_type == cls):
                if isinstance(constructor.boost_python_export, bool):
                    constructor.boost_python_export = False
    # for cls in asg['class ::std::set'].specializations(partial = False):
    #     for method in cls.methods():
    #         if method.localname in ['swap', 'key_comp', 'value_comp', 'get_allocator']:
    #             if isinstance(method.boost_python_export, bool):
    #                 method.boost_python_export = False
    #     for constructor in cls.constructors():
    #         if not(constructor.nb_parameters == 0 or constructor.nb_parameters == 1 and constructor.parameters[0].qualified_type.unqualified_type == cls):
    #             if isinstance(constructor.boost_python_export, bool):
    #                 constructor.boost_python_export = False
    # for cls in asg['class ::std::less'].specializations(partial = False):
    #     cls.boost_python_export = False
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
    # for spec in asg['class ::std::basic_ostream'].specializations(partial = False):
    #     spec.boost_python_export = False
    for cls in ['class ::std::basic_ostream< char, struct ::std::char_traits< char > >',
                'class ::std::basic_streambuf< char, struct ::std::char_traits< char > >',
                'class ::std::basic_filebuf< char, struct ::std::char_traits< char > >',
                'class ::std::basic_ifstream< char, struct ::std::char_traits< char > >']:
        asg[cls].boost_python_export = False
    return asg

if __name__ == "__main__":
    asg = autowig.AbstractSemanticGraph()
    prefix = Path('../../build-scons/include/stat_tool')
    autowig.parser.plugin = 'pyclanglite'
    asg = autowig.parser(asg, prefix.walkfiles('*.h*'),
                              ['-x', 'c++', '-std=c++11',  '-I' + str(prefix.parent.abspath()),
                               '-I' + sys.prefix + '/include', '-I' + sys.prefix + '/include/python2.7'],
                              silent = True)
    asg = controller(asg)
    autowig.generator.plugin = 'boost_python_internal'
    wrappers = autowig.generator(asg,
                             module = '../py/wrapper/_stat_tool.cpp',
                             decorator = '../py/structure_analysis/stat_tool/_stat_tool.py',
                             prefix = 'wrapper_')
    wrappers.write()