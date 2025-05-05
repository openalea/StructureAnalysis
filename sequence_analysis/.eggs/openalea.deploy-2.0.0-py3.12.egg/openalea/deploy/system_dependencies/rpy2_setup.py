
import os, os.path, sys, shutil, re, itertools
from collections import namedtuple
from distutils.command.build_ext import build_ext as _build_ext
from distutils.command.build import build as _build

from distutils.core import setup
from distutils.core import Extension

if sys.platform == "win32":
    def GetShortPathName(path):
        import ctypes
        nsize = 2048
        out = ctypes.create_unicode_buffer(nsize)
        short = None
        if ctypes.windll.kernel32.GetShortPathNameW(path, ctypes.pointer(out), nsize):
            return ctypes.wstring_at(path)
        else:
            print ctypes.windll.kernel32.GetLastError()


pack_name = 'rpy2'
pack_version = __import__('rpy').__version__

package_prefix='.'
if sys.version_info >= (3,):
    print("Using 2to3 to translate Python2-only idioms into Python3 code. Please wait...")
    # Python 3 and we need to translate code
    package_prefix = os.path.join('build', 'python3_rpy')
    from distutils import filelist, dir_util, file_util, util#, log
    #log.set_verbosity(1)
    fl = filelist.FileList()
    tmp = open("MANIFEST.in")
    for line in tmp:
        line = line.rstrip()
        if line != '':
            fl.process_template_line(line)
    tmp.close()
    dir_util.create_tree(package_prefix, fl.files)
    outfiles_2to3 = []
    #dist_script = os.path.join("build", "src", "distribute_setup.py")
    for f in fl.files:
        outf, copied = file_util.copy_file(f, os.path.join(package_prefix, f), 
                                           update=1)
        if copied and outf.endswith(".py"): #and outf != dist_script:
            outfiles_2to3.append(outf)
        if copied and outf.endswith('api_tests.txt'):
            # XXX support this in distutils as well
            from lib2to3.main import main
            main('lib2to3.fixes', ['-wd', os.path.join(package_prefix, 
                                                       'tests', 'api_tests.txt')])

    util.run_2to3(outfiles_2to3)

    # arrange setup to use the copy
    sys.path.insert(0, package_prefix)
    src_root = package_prefix
    print('done.')
else:
    from distutils.core import setup    
from distutils.core import Extension


class build(_build):
    user_options = _build.user_options + \
        [
        #('r-autoconfig', None,
        # "guess all configuration paths from " +\
        #     "the R executable found in the PATH " +\
        #     "(this overrides r-home)"),
        ('r-home=', None, 
         "full path for the R home to compile against " +\
             "(see r-autoconfig for an automatic configuration)"),
        ('r-home-lib=', None,
         "full path for the R shared lib/ directory " +\
             "(<r-home>/lib otherwise)"),
        ('r-home-modules=', None,
         "full path for the R shared modules/ directory " +\
             "(<r-home>/modules otherwise)") 
        ]
    boolean_options = _build.boolean_options #+ \
        #['r-autoconfig', ]


    def initialize_options(self):
        _build.initialize_options(self)
        self.r_autoconfig = None
        self.r_home = None
        self.r_home_lib = None
        self.r_home_modules = None

class build_ext(_build_ext):
    """
    -DRPY_STRNDUP          : definition of strndup()
    -DRPY_VERBOSE
    -DRPY_DEBUG_PRESERV
    -DRPY_DEBUG_PROMISE    : evaluation of promises
    -DRPY_DEBUG_OBJECTINIT : initialization of PySexpObject
    -DRPY_DEBUG_CONSOLE    : console I/O
    -DRPY_DEBUG_COBJECT    : SexpObject passed as a CObject
    -DRPY_DEBUG_GRDEV
    """
    user_options = _build_ext.user_options + \
        [
        #('r-autoconfig', None,
        #  "guess all configuration paths from " +\
        #      "the R executable found in the PATH " +\
        #      "(this overrides r-home)"),
        ('r-home=', None, 
         "full path for the R home to compile against " +\
             "(see r-autoconfig for an automatic configuration)"),
        ('r-home-lib=', None,
         "full path for the R shared lib/ directory" +\
             "(<r-home>/lib otherwise)"),
        ('r-home-modules=', None,
         "full path for the R shared modules/ directory" +\
             "(<r-home>/modules otherwise)")]

    boolean_options = _build_ext.boolean_options #+ \
        #['r-autoconfig', ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.r_autoconfig = None
        self.r_home = None
        self.r_home_lib = None
        self.r_home_modules = None

    def finalize_options(self):
        self.set_undefined_options('build',
                                   #('r_autoconfig', 'r_autoconfig'),
                                   ('r_home', 'r_home'))
        _build_ext.finalize_options(self) 
        if self.r_home is None:
            tmp = os.popen("R RHOME")
            self.r_home = tmp.readlines()
            tmp.close()
            if len(self.r_home) == 0:
                raise SystemExit("Error: Tried to guess R's HOME but no R command in the PATH.")

    #Twist if 'R RHOME' spits out a warning
            if self.r_home[0].startswith("WARNING"):
                self.r_home = self.r_home[1]
            else:
                self.r_home = self.r_home[0]
            #self.r_home = [self.r_home, ]

        if self.r_home is None:
            raise SystemExit("Error: --r-home not specified.")
        else:
            self.r_home = self.r_home.split(os.pathsep)

        rversions = []
        for r_home in self.r_home:
            r_home = r_home.strip()
        rversion = get_rversion(r_home)
        if cmp_version(rversion[:2], [2, 8]) == -1:
            raise SystemExit("Error: R >= 2.8 required.")
        rversions.append(rversion)

        config = RConfig()
        for about in ('--ldflags', '--cppflags', 
                      'LAPACK_LIBS', 'BLAS_LIBS'):
            config += get_rconfig(r_home, about)            

        print(config.__repr__())

        self.include_dirs.extend(config._include_dirs)
        self.libraries.extend(config._libraries)
        self.library_dirs.extend(config._library_dirs)        

        if self.r_home_modules is None:
            self.library_dirs.extend([os.path.join(r_home, 'modules'), ])
        else:
            self.library_dirs.extends([self.r_home_modules, ])

        # for e in self.extensions:
        #    self.extra_link_args.extra_link_args(config.extra_link_args)
        #    e.extra_compile_args.extend(extra_compile_args)

    def run(self):
        _build_ext.run(self)



def get_rversion(r_home):
    r_exec = os.path.join(r_home, 'bin', 'R')
    # Twist if Win32
    if sys.platform == "win32":
        rp = os.popen3('"'+r_exec+'" --version')[2]
    else:
        rp = os.popen('"'+r_exec+'" --version')
    rversion = rp.readline()
    #Twist if 'R RHOME' spits out a warning
    if rversion.startswith("WARNING"):
        rversion = rp.readline()
    m = re.match('^R version ([^ ]+) .+$', rversion)
    rversion = m.groups()[0]
    rversion = rversion.split('.')
    rversion[0] = int(rversion[0])
    rversion[1] = int(rversion[1])
    rp.close()
    return rversion

def cmp_version(x, y):
    if (x[0] < y[0]):
        return -1
    if (x[0] > y[0]):
        return 1
    if (x[0] == y[0]):
        if len(x) == 1 or len(y) == 1:
            return 0
        return cmp_version(x[1:], y[1:])

class RConfig(object):
    _include_dirs = None
    _libraries = None
    _library_dirs = None 
    _extra_link_args = None
    _frameworks = None
    _framework_dirs = None

    def __init__(self,
                 include_dirs = tuple(), libraries = tuple(),
                 library_dirs = tuple(), extra_link_args = tuple(),
                 frameworks = tuple(),
                 framework_dirs = tuple()):
        for k in ('include_dirs', 'libraries', 
                  'library_dirs', 'extra_link_args'):
            v = locals()[k]
            if not isinstance(v, tuple):
                if isinstance(v, str):
                    v = [v, ]
            v = tuple(set(v))
            self.__dict__['_'+k] = v
        # frameworks are specific to OSX
        for k in ('framework_dirs', 'frameworks'):
            v = locals()[k]
            if not isinstance(v, tuple):
                if isinstance(v, str):
                    v = [v, ]
            v = tuple(set(v))
            self.__dict__['_'+k] = v
            self.__dict__['_'+'extra_link_args'] = tuple(set(v + self.__dict__['_'+'extra_link_args']))

    @staticmethod
    def from_string(string, allow_empty = False):
        possible_patterns = ('^-L(?P<library_dirs>[^ ]+)$',
                             '^-l(?P<libraries>[^ ]+)$',
                             '^-I(?P<include_dirs>[^ ]+)$',
                             '^(?P<framework_dirs>-F[^ ]+?)$',
                             '^(?P<frameworks>-framework [^ ]+)$',
                             '^(?P<extra_link_args>-Wl[^ ]+)$')
        pp = [re.compile(x) for x in possible_patterns]
        # sanity check of what is returned into rconfig
        rconfig_m = None        
        span = (0, 0)
        rc = RConfig()
        
        for substring in re.split('(?<!-framework) ', string):
            ok = False
            if not substring:
                continue
            for pattern in pp:
                rconfig_m = pattern.match(substring)
                if rconfig_m is not None:
                    rc += RConfig(**rconfig_m.groupdict())
                    span = rconfig_m.span()
                    ok = True
                    break
                elif rconfig_m is None:
                    if allow_empty and (rconfig == ''):
                        print(cmd + '\nreturned an empty string.\n')
                        rc += RConfig()
                        ok = True
                        break
                    else:
                        # if the configuration points to an existing library, 
                        # use it
                        if os.path.exists(string):
                            rc += RConfig(library = substring)
                            ok = True
                            break
            if not ok:
                raise ValueError('Invalid substring\n' + substring 
                                 + '\nin string\n' + string)
        return rc
            
    def __repr__(self):
        s = 'Configuration for R as a library:' + os.linesep
        s += os.linesep.join(
            ['  ' + x + ': ' + self.__dict__['_'+x].__repr__() \
                 for x in ('include_dirs', 'libraries',
                           'library_dirs', 'extra_link_args')])
        s += os.linesep + ' # OSX-specific (included in extra_link_args)' + os.linesep 
        s += os.linesep.join(
            ['  ' + x + ': ' + self.__dict__['_'+x].__repr__() \
                 for x in ('framework_dirs', 'frameworks')]
            )
        
        return s

    def __add__(self, config):
        assert isinstance(config, RConfig)
        res = RConfig(include_dirs = self._include_dirs + \
                          config._include_dirs,
                      libraries = self._libraries + config._libraries,
                      library_dirs = self._library_dirs + \
                          config._library_dirs,
                      extra_link_args = self._extra_link_args + \
                          config._extra_link_args,
                      frameworks = self._frameworks + config._frameworks,
                      framework_dirs = self._framework_dirs + config._framework_dirs)
        return res


def get_rconfig(r_home, about, allow_empty = False):
    if sys.platform == "win32":
        # parse the make conf file      
        arch = 'i386'
        inc_dirs     = os.path.join(r_home, 'include')
        arch_inc_dir = os.path.join(r_home, 'etc', arch)
        print os.path.join(r_home, 'include')
        print arch_inc_dir
        d = get_makeconf_dict(r_home, arch)
        if about == "--cppflags":
            rconfig = d.get("ALL_CPPFLAGS") + r' -I"%s" -I"%s"'%(inc_dirs, arch_inc_dir)
        elif about == "--ldflags":
            rconfig = d.get("ALL_LIBS") + " -lreadline"
        elif about == "LAPACK_LIBS":
            rconfig = d.get("LAPACK_LIBS")        
        elif about == "BLAS_LIBS":
            rconfig = d.get("BLAS_LIBS")
    else:
        r_exec = os.path.join(r_home, 'bin', 'R')
        cmd = '"'+r_exec+'" CMD config '+about
        rp = os.popen(cmd)
        rconfig = rp.readline()
        #Twist if 'R RHOME' spits out a warning
        if rconfig.startswith("WARNING"):
            rconfig = rp.readline()
        rp.close()
    print rconfig
    rconfig = rconfig.strip()
    rc = RConfig.from_string(rconfig, allow_empty = allow_empty)
    return rc

    
####################
# PATCHED SECTION  #
####################   
import re
from collections import OrderedDict
ass_re = re.compile(r"^\s*([\-\w]+)\s*=\s*(.*)$")
app_re = re.compile(r"^\s*(\[\-\w]+)\s*\+=\s*(.*)$")
var_re = re.compile(r"\$\(([\-\w]*)\)$")
var_re2 = re.compile(r"(.*)\$\(([\-\w]*)\)(.*)")
def get_makeconf_dict(r_home, arch):
    makeconf = os.path.join(r_home, 'etc', arch, 'Makeconf')
    return parse_make(makeconf, initial_vars={"R_HOME": r_home} )    

def parse_make(makefile, initial_vars=None):	
    if isinstance(makefile, file):
        makefile.seek(0)
        txt = makefile.readlines()
    else:
        with open(makefile) as f:
            txt = f.readlines()
    
    # -- gather variables --
    vardict = {}
    # vardeps = [] #: variables interdependencies    
    
    for l in txt:
        if "=" in l and not l.startswith("#"):

            # ignore commented equal signs
            if "#" in l and l.index("#") < l.index("="): continue
            
			# we do not do rule expansion so $* and relatives are not recognized
            l = l.replace("$*","")
			
            assignmt = ass_re.match(l)
            appendmt = app_re.match(l)
            if assignmt:                
                k = assignmt.group(1)
                v = assignmt.group(2)

                #print k,v
                k = k.strip()
                v = v.strip()          
                vardict[k] = v
                
                se = var_re.findall(v)
                # if se:
                    # vardeps.append( (k, se) )
            elif(appendmt):
                k = appendmt.group(1)
                v = appendmt.group(2)
                if k in vardict:
                    vardict[k] += " "+v
                else:
                    vardict[k] = v
                    
    if initial_vars:
        vardict.update(initial_vars)
            
    # expand variables
    for var, value in vardict.copy().iteritems():
        expand = var_re2.search(value) #something to expand
        while (expand):
            value  = expand.group(1) + vardict.get(expand.group(2), "") + expand.group(3)  
            expand = var_re2.search(value) #something to expand
        vardict[var] = value

    return vardict    
#####################
# /PATCHED SECTION  #
##################### 

def getRinterface_ext():
    #r_libs = [os.path.join(RHOME, 'lib'), os.path.join(RHOME, 'modules')]
    r_libs = []
    extra_link_args = []

    #FIXME: crude way (will break in many cases)
    #check how to get how to have a configure step
    define_macros = []

    if sys.platform == 'win32':
        define_macros.append(('Win32', 1))
        define_macros.append(('NOGDI', None))               #avoids redefining ERROR
        define_macros.append(('_RPY_STRNDUP_', None))       #enables own definition of strndup
        define_macros.append(('_UINTPTR_T_DEFINED', None))  #avoids definition of uintptr_t by stdint.h
        define_macros.append(('uintptr_t', "unsigned int")) #avoids redefinition of uintptr_t by Rinterface.h
    else:
        define_macros.append(('R_INTERFACE_PTRS', 1))
        define_macros.append(('HAVE_POSIX_SIGJMP', 1))
        define_macros.append(('RIF_HAS_RSIGHAND', 1))
        
    define_macros.append(('CSTACK_DEFNS', 1))


    if sys.byteorder == 'big':
        define_macros.append(('RPY_BIGENDIAN', 1))
    else:
        pass


    include_dirs = []
    
    rinterface_ext = Extension(
            name = pack_name + '.rinterface._rinterface',
            sources = [ \
            #os.path.join('rpy', 'rinterface', 'embeddedr.c'), 
            #os.path.join('rpy', 'rinterface', 'r_utils.c'),
            #os.path.join('rpy', 'rinterface', 'buffer.c'),
            #os.path.join('rpy', 'rinterface', 'sequence.c'),
            #os.path.join('rpy', 'rinterface', 'sexp.c'),
            os.path.join(package_prefix,
                         'rpy', 'rinterface', '_rinterface.c')
                       ],
            depends = [os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'embeddedr.h'), 
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'r_utils.h'),
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'buffer.h'),
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'sequence.h'),
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'sexp.h'),
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'rpy_rinterface.h'),
                       os.path.join(package_prefix,
                                    'rpy', 'rinterface', 'rpy_device.h')
                       ],
            include_dirs = [os.path.join(package_prefix,
                                         'rpy', 'rinterface'),] + include_dirs,
            libraries = ['R'],
            library_dirs = r_libs,
            define_macros = define_macros,
            runtime_library_dirs = r_libs,
            extra_compile_args=[],
            #extra_link_args = extra_link_args
            )
            
    rpy_device_ext = Extension(
        pack_name + '.rinterface.rpy_device',
            [
            os.path.join('rpy', 'rinterface', 'rpy_device.c'),
             ],
            include_dirs = include_dirs + 
                            [os.path.join('rpy', 'rinterface'), ],
            libraries = ['R', ],
            library_dirs = r_libs,
            define_macros = define_macros,
            runtime_library_dirs = r_libs,
            #extra_compile_args=['-O0', '-g'],
            extra_link_args = extra_link_args
        )

    return [rinterface_ext, rpy_device_ext]

if __name__ == '__main__':
    rinterface_exts = []
    ri_ext = getRinterface_ext()
    rinterface_exts.append(ri_ext)

    pack_dir = {pack_name: os.path.join(package_prefix, 'rpy')}

    import distutils.command.install
    for scheme in distutils.command.install.INSTALL_SCHEMES.values():
        scheme['data'] = scheme['purelib']

    setup(
        #install_requires=['distribute'],
        cmdclass = {'build': build,
                    'build_ext': build_ext},
        name = pack_name,
        version = pack_version,
        description = "Python interface to the R language",
        url = "http://rpy.sourceforge.net",
        license = "AGPLv3.0 (except rpy2.rinterface: LGPL)",
        author = "Laurent Gautier",
        author_email = "lgautier@gmail.com",
        ext_modules = rinterface_exts[0],
        package_dir = pack_dir,
        packages = [pack_name,
                    pack_name + '.rlike',
                    pack_name + '.rlike.tests',
                    pack_name + '.rinterface',
                    pack_name + '.rinterface.tests',
                    pack_name + '.robjects',
                    pack_name + '.robjects.tests',
                    pack_name + '.robjects.lib',
                    pack_name + '.interactive',
                    pack_name + '.interactive.tests'
                    ],
        classifiers = ['Programming Language :: Python',
                       'Programming Language :: Python :: 2',
                       'Programming Language :: Python :: 3',
                       'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
                       'License :: OSI Approved :: GNU Affero General Public License v3',
                       'Intended Audience :: Developers',
                       'Intended Audience :: Science/Research',
                       'Development Status :: 5 - Production/Stable'
                       ],
        data_files = [(os.path.join('rpy2', 'images'), 
                       [os.path.join('doc', 'source', 'rpy2_logo.png')])],

        #[pack_name + '.rinterface_' + x for x in rinterface_rversions] + \
            #[pack_name + '.rinterface_' + x + '.tests' for x in rinterface_rversions]
        )

