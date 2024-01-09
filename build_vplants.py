#-*-python-*-

import os, sys

#cmd = "python setup.py install"
cmd = "scons install"
if os.name == 'nt':
    compiler = 'mingw'
    cmd = 'python setup.py build_ext --scons-ext-param="compiler=%s" install'%(compiler,)  
    print cmd
else:
    cmd = "python setup.py install --install-lib='/home/bolivier/python-lib' --install-dyn-lib='/home/bolivier/python-lib/lib/'"

# dependencies dictionnary

deps = {"tool": [],
        #"plantgl" : [],
        "stat_tool": ["tool"],
        "sequence_analysis": ["stat_tool"],
        "amlobj" :['tool'],
        #"mtg":["plantgl", "amlobj"],
        "mtg":["amlobj"],
        "tree":[],
        "tree_matching":["mtg", "stat_tool"],
        #"msvoxel":["mtg"],
        #"flux":["mtg"],
        "aml":["sequence_analysis", "tree_matching"], #, "flux", "msvoxel"],
        "tree_statistic": ["tree", "sequence_analysis"],
        "all" : ["tree_statistic", "fractalysis", "aml"],
        "all" : ["fractalysis", "aml"],
        }

dirs = {"plantgl" : "PlantGL",}




def compute_deps(target, deps):
    """ Return a list of dependencies of target in deps """

    l = [target]
    if(deps.has_key(target)):
        for d in deps[target]:
            l += compute_deps(d, deps)
            
    return l

    

# current directpry
cdir = os.path.abspath(os.curdir)

# get target

if __name__ == "__main__":
    
    finished = set() # target already executed
    targets = sys.argv[1:]
    if(not targets) : targets = ["all"]

    for t in targets:
        depends = compute_deps(t, deps)
        depends.reverse()

        for d in depends:
            if(d in finished) : continue
            if d == 'all': continue

            print "Installing :", d

            directory = dirs.get(d, d)
            os.chdir(os.path.join(cdir, directory))
            status = os.system(cmd)
            if(status > 0):
                print "Error in building or Installing ",d
                print "Please, add specific options in a file options.py"
                sys.exit(1)
            finished.add(d)

        
