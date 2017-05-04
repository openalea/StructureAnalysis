# rm src/py/wrapper/*.cpp
# rm src/py/wrapper/*.h
# source activate vplants && scons; source deactivate vplants
# cd src/cpp && source activate statiskit-dev && python AutoWIG.py; source deactivate statiskit-dev cd ../..
# scons

import autowig
from path import Path
import sys

def controller(asg):
    # for spec in asg['class ::std::basic_ostream'].specializations(partial = False):
    #     spec.boost_python_export = False
    for cls in ['class ::std::basic_ostream< char, struct ::std::char_traits< char > >',
                'class ::std::basic_streambuf< char, struct ::std::char_traits< char > >',
                'class ::std::basic_filebuf< char, struct ::std::char_traits< char > >',
                'class ::std::basic_ifstream< char, struct ::std::char_traits< char > >']:
        asg[cls].boost_python_export = False
    return autowig.controller(asg)    

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