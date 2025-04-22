# -*-python-*-

from openalea.sconsx import config, environ
import os, platform


pj = os.path.join
ALEASolution = config.ALEASolution

options = Variables(['../options.py', 'options.py'], ARGUMENTS)
# Firstly get options in ../options.py and then in options.py and finally in ARGUMENTS
options.Add(BoolVariable('with_efence', 'build with efence library', False))
options.Add(BoolVariable('with_test', 'build with efence library', False))


tools = ['boost_python', 'openalea.stattool']

env = ALEASolution(options, tools)
env.AppendUnique(CXXFLAGS=['-x', 'c++', '-std=c++14'])

if env['with_efence']:
    env.AppendUnique(LIBS=['efence'])

if (platform.system() != 'Windows' and
    os.environ.get('CC') and
    os.environ.get('CXX')):

    env.AppendUnique(CFLAGS=["-std=c14"])
    if (platform.system() == 'Darwin'):
        env.AppendUnique(CXXFLAGS=['-stdlib=libc++'])
    conda_prefix = os.environ['CONDA_PREFIX']
    conda_bin = pj(os.environ['CONDA_PREFIX'], 'bin')
    
    # To work with conda toolchain
    env['AR'] = os.environ['AR']
    env['AS'] = os.environ['AS']
    env['CC'] = pj(conda_bin, os.environ['CC'])
    env['CXX'] = pj(conda_bin, os.environ['CXX'])
    env.PrependUnique(
        CPPPATH=['%s/include'%(conda_prefix)],
        CCFLAGS=['-fvisibility=hidden'],
        LIBPATH=['%s/lib'%(conda_prefix)])

# Build stage
prefix = env['build_prefix']
seqlib = SConscript(pj(prefix,"src/cpp/SConscript"),
        exports='env')
SConscript(pj(prefix,"src/wrapper/SConscript"),
        exports='env')

if bool(env['with_test']):
    SConscript(pj(prefix,"test/cpp/SConscript"), exports="env seqlib")

Default("build")

