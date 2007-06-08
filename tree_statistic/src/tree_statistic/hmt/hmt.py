"""Hidden Markov tree models
"""
import string
import stat_tool, ctree, ctrees, chmt, trees

VariableType=stat_tool.VariableType
FormatError=stat_tool.FormatError
CharacteristicType=ctree.Characteristic

class HiddenMarkovTree:
    """An implementation of the hidden Markov out-trees with conditionally 
    independent children states given their parent."""

    def __init__(self, arg, aliasing=False):
        """Initialize a Hmt by copy or by reading into a file.

            Usage:  H=HiddenMarkovTree("file_name")
                    H=HiddenMarkovTree(HiddenMarkovTree)."""
        ## Aliasing is used to make an alias between argument and
        ## self.__chmt -> useful for the connection between HiddenMarkovTree
        ## and HiddenMarkovTreeData
        if type(arg)==str:
            # arg is expectedly a file name...
            try:
                self.__chmt=chmt.HmtAsciiRead(arg)
            except RuntimeError, error:
                raise FormatError, error
            nbvariables=self.__chmt.NbInt()+self.__chmt.NbFloat()
            self._attributes=["Variable " + str(i) for i in range(nbvariables)]
        elif issubclass(arg.__class__, HiddenMarkovTree):
            # ... or a HiddenMarkovTree object...
            if aliasing:
                self.__chmt=arg.__chmt
            else:
                self.__chmt=chmt.CiHmot(arg.__chmt)
            self._attributes=list(arg._attributes)
        elif issubclass(arg.__class__, chmt.CiHmot):
            # ... or a CiHmot object
            if aliasing:
                self.__chmt=arg
            else:
                self.__chmt=chmt.CiHmot(arg)
            nbvariables=self.__chmt.NbInt()+self.__chmt.NbFloat()
            self._attributes=["Variable " + str(i) for i in range(nbvariables)]
        else:
            msg="bad argument type: "+str(type(arg))
            raise TypeError, msg

    def Display(self, ViewPoint=None, Detail=None, TreeId=None, 
                NbStateTrees=2, StateTrees="GeneralizedViterbi"):
        """Display HiddenMarkovTree object.
        
        Usage: Display()
               Display(ViewPoint="Data", Detail=2)
               Display(ViewPoint="StateProfile", TreeId=0)"""
        if ViewPoint is None:
            wp="Data"
        else:
            wp=ViewPoint
        if type(wp)!=str:        
            msg='bad type for argument "ViewPoint": '+str(type(wp))
            raise TypeError, msg
        if (wp.upper()=="DATA"):
            if (TreeId is None):
                if ((Detail==1) or (Detail is None)):
                    exhaustive=False
                elif Detail==2:
                    exhaustive=True
                elif type(Detail)!=int:
                    msg="Bad type for 'Detail' argument:"+str(type(Detail)) \
                    +" - expecting type 'int'"
                    raise TypeError, msg
                else:
                    msg="Bad value for 'Detail' argument:"+str(Detail) \
                    +" - expecting 1 or 2"
                    raise ValueError, msg
                print self.__chmt.Display(exhaustive)
            else:
                msg="Argument 'TreeId' is not required when using "+\
                'Display() or Display(ViewPoint="Data")'
                raise ValueError, msg
        elif (wp.upper()=="STATEPROFILE"):
            if (Detail is None):
                if (TreeId is None):
                    raise TypeError, "Argument 'TreeId' is mandatory in "+\
                            'Display(ViewPoint="StateProfile", TreeId)'
                if (type(NbStateTrees)!=int):
                    msg='bad type for argument "NbStateTrees": '\
                        +str(type(NbStateTrees))
                    raise TypeError, msg
                try:
                    ## chmt_data=self.__chmt.ExtractData()
                    res=self.__chmt.StateProfile(0, NbStateTrees, TreeId)
                except RuntimeError, error:
                    raise FormatError, error
                else:
                    profile_nature=0
                    for o in res:
                        if (type(o)==str):
                            print o
                        else:
                            if profile_nature==0:
                                # smoothing
                                nb_int=self.__chmt.NbInt()
                                nb_float=self.__chmt.NbFloat()
                                nbvariables=o.NbInt()+\
                                    o.NbFloat()-self.NbStates()-1
                                attributes= ["Optimal State"]+\
                                    self._attributes[:nb_int+1]+\
                                    ["State "+str(v) for v in range(self.NbStates())]+\
                                    ["CalEntropy", "MalEntropy", "PalEntropy"]+\
                                    self._attributes[nb_int+1:]
                                smoothed=HiddenMarkovTreeData(o,
                                                              markov=self.__chmt, 
                                                              aliasing=True, 
                                                              attribute_names=attributes)
                                smoothed.Tree(0).Display()
                            elif profile_nature==1:
                                # viterbi upward downward
                                nbvariables=o.NbInt()+\
                                    o.NbFloat()-self.NbStates()-1
                                attributes= ["Optimal State"]+\
                                    ["Variable "+str(v) for v in range(nbvariables)]+\
                                    ["State "+str(v) for v in range(self.NbStates())]
                                viterbi_ratios=\
                                    HiddenMarkovTreeData(o, markov=self.__chmt, 
                                                        aliasing=True, 
                                                        attribute_names=attributes)
                                viterbi_ratios.Tree(0).Display()
                            elif profile_nature==2:
                                # suboptimal state trees
                                attributes=\
                                    ["Rank "+str(v+1) for v in range(o.NbInt())]
                                suboptimal_states=\
                                    HiddenMarkovTreeData(o, markov=self.__chmt, 
                                                        aliasing=True, 
                                                        attribute_names=attributes)
                                sst=suboptimal_states.Tree(0)
                                sst._set_int_value_type(0)
                                sst.Display()
                            else:    
                                chmt=HiddenMarkovTreeData(o, markov=self.__chmt, 
                                                          aliasing=True)
                                chmt.Tree(0).Display()
                                msg="Incorrect type of state profile: " + \
                                   str(profile_nature)
                            profile_nature+=1
            else:
                msg="Argument 'Detail' is not required when using "+\
                'Display(ViewPoint="StateProfile")'
                raise ValueError, msg
        else:
            msg='bad value for argument "ViewPoint": '+ViewPoint
            raise ValueError, msg
        
        
    def ExtractData(self):
        """Extract the 'data' part of the HiddenMarkovTree."""
        try:
            chmt_data=self.__chmt.ExtractData()
        except RuntimeError, error:
            raise FormatError, error
        else:
            chmt_data=chmt_data.StateTrees()
            return HiddenMarkovTreeData(chmt_data, markov=None, aliasing=True)

    def ExtractDistribution(self, type, arg1, arg2=None):
        """Extract a distribution from the HiddenMarkovTree.
        
        Usage:  ExtractDistribution("Observation", variable, state)"""
        try:
            cdistribution_data=self.__chmt.ExtractDistribution(arg1, arg2)
        except RuntimeError, error:
            raise FormatError, error
        else:
            return Distribution(cdistribution_data)

    def ExtractPlot(self, TreeId, ViewPoint):
        """Extract a tree with state or entropy profile.
        
        Usage:  ExtractPlot(TreeId, ViewPoint="StateProfile")"""
        if (TreeId is None):
            raise TypeError, "Argument 'TreeId' is mandatory in "+\
                    'ExtractPlot(TreeId, ViewPoint="StateProfile")'
        try:
            res=self.__chmt.StateProfile(0, 2, TreeId)
        except RuntimeError, error:
            raise FormatError, error
        else:
            # string result
            resprint=[]
            profile_nature=0
            for o in res:
                if (type(o)!=str):
                    if profile_nature==0:
                        # smoothing
                        nb_int=self.__chmt.NbInt()
                        nb_float=self.__chmt.NbFloat()
                        nbvariables=o.NbInt()+\
                            o.NbFloat()-self.NbStates()-1
                        attributes= ["Optimal State"]+\
                            self._attributes[:nb_int+1]+\
                            ["State "+str(v) for v in range(self.NbStates())]+\
                            ["CalEntropy", "MalEntropy", "PalEntropy"]+\
                            self._attributes[nb_int+1:]
                        smoothed=HiddenMarkovTreeData(o, 
                                                      markov=self.__chmt, 
                                                      aliasing=True, 
                                                      attribute_names=attributes)
                        # smoothed.Tree(0).Display()
                        profile_nature+=1
                        return smoothed
##                    elif profile_nature==1:
##                        # maximum state probabilities = Viterbi upward-downward
##                        nb_int=self.__chmt.NbInt()
##                        nb_float=self.__chmt.NbFloat()
##                        nbvariables=o.NbInt()+\
##                            o.NbFloat()-self.NbStates()-1
##                        attributes= ["Optimal State"]+\
##                            self._attributes[:nb_int+1]+\
##                            ["State "+str(v) for v in range(self.NbStates())]+\
##                            self._attributes[nb_int+1:]
##                        viterbiud=HiddenMarkovTreeData(o, 
##                                                       markov=self.__chmt, 
##                                                       aliasing=True, 
##                                                       attribute_names=attributes)
##                        viterbiud.Tree(0).Display()
##                        profile_nature+=1
##                else:
##                    print o                    
        
    def NbStates(self):
        """Return the number of hidden states of the HiddenMarkovTree."""
        return self.__chmt.NbStates()
    
    def Plot(self, ViewPoint="Observation", variable=0, Title="", 
             TreeId=0, EndVertex=None):
        """Graphical output using Gnuplot.py.        
        
        Usage:  Plot(ViewPoint="StateProfile", TreeId=0, EndVertex=10)
                Plot(ViewPoint="FirstOccurrenceRoot", variable=0)
        Other possible values for ViewPoint: 
            FirstOccurrenceLeaves
            SojournSize
            Counting
            Observation
            """
        if type(ViewPoint)!=str:
            msg='bad type for argument "ViewPoint": '\
              +"type 'str' expected"
            raise TypeError, msg
        import os
        is_observation=False
        if ViewPoint.upper()=="FIRSTOCCURRENCEROOT":
            ftype=CharacteristicType.FIRST_OCCURRENCE_ROOT
        elif ViewPoint.upper()=="FIRSTOCCURRENCELEAVES":
            ftype=CharacteristicType.FIRST_OCCURRENCE_LEAVES
        elif ViewPoint.upper()=="SOJOURNSIZE":
            ftype=CharacteristicType.SOJOURN_SIZE
        elif ViewPoint.upper()=="COUNTING":
            ftype=CharacteristicType.NB_ZONES
        elif ViewPoint.upper()=="OBSERVATION":
            ftype=0 #CharacteristicType.OBSERVATION
            nb_windows=self.__chmt.NbStates()
            is_observation=True
        elif ViewPoint.upper()!="STATEPROFILE":
            msg='bad value for argument "ViewPoint": '+ViewPoint
            raise ValueError, msg
        if ViewPoint.upper()!="STATEPROFILE":
            try:
                w=self.__chmt.NbValues(variable)
            except RuntimeError:
                msg="variable index out of range: "+str(variable)
                raise IndexError, msg
            if not(is_observation):
                nb_windows=w
                
            # part of the filename which identifies the involved variable
            ref_file_id=str(variable+1)
            # part of the filename which identifies the graph to be displayed
            file_id=ref_file_id
            if not(is_observation):
                file_id+=str(ftype+2)
            # in the special case of observation distributions, 
            # the suffix can be "0" or void, depending on the number of values
            prefix="ftmp"
            file_created=False
            file_list=[]
            # find a non existing file name
            while not file_created:
                try:
                    cfile=open(prefix+'01.dat','r')
                except IOError:
                    file_created=True
                else:
                    import random
                    prefix+=str(random.randint(1,9))
            try:
                self.__chmt.Plot(os.getcwd()+os.sep+prefix, Title)
                # build the list of the created files: 
                for var in range(self.__chmt.NbInt()+1):
                    for char in [str(c) for c in range(5)]+[""]:
                        filename=prefix+str(var)+char
                        try:
                            tmpfile=open(filename+'.plot', 'r')
                        except IOError:
                            pass
                        else:
                            tmpfile.close()
                            # add the .plot and .print files
                            file_list+=[filename+extension
                                for extension in [".plot", ".print"]]
                            if (char=='2'):
                                # add the .dat file
                                file_list+=[prefix+str(var)+"1.dat"]
                            elif  (char==''):
                                file_list+=[prefix+str(var)+".dat"]
                            # get the correct file_id in the case 
                            # of observation distributions
                            if (((var==variable+1) and is_observation) 
                                and (char=='0')):
                                file_id+="0"
            except RuntimeError, f:
                for tmpfile in file_list:
                   os.remove(tmpfile)
                raise FormatError, f
        else: # ViewPoint.upper()=="STATEPROFILE"
            nb_windows=4
            if (type(TreeId)!=int):
                raise TypeError, "bad tree index type: "+str(type(TreeId))
            elif (TreeId < 0):
                raise IndexError, "bad tree index value: "+str(variable)
            elif (type(EndVertex)!=int):
                msg="bad type for vertex identifier: " + \
                         str(type(EndVertex))
                raise TypeError, msg
            ref_file_id=""
            file_id=ref_file_id
            prefix="ftmp"
            file_created=False
            file_list=[]
            # find a non existing file name
            while not file_created:
                try:
                    cfile0=open(prefix+'0.dat','r')
                except IOError:
                    try:
                        cfile1=open(prefix+'1.dat','r')
                    except IOError:
                        file_created=True
                    else:
                        import random
                        prefix+=str(random.randint(1,9))
                else:
                    import random
                    prefix+=str(random.randint(1,9))
            try:
                self.__chmt.StateProfilePlot(os.getcwd()+os.sep+prefix, 
                                             Title,
                                             TreeId,
                                             EndVertex)
                # build the list of the created files: 
                for char in [str(c) for c in range(2)]+[""]:
                    filename=prefix+char
                    try:
                        tmpfile=open(filename+'.dat', 'r')
                    except IOError:
                        pass
                    else:
                        tmpfile.close()
                        file_list+=[filename+".dat"]
                file_list+=[prefix+extension
                    for extension in [".plot", ".print"]]
            except RuntimeError, f:
                for tmpfile in file_list:
                   os.remove(tmpfile)
                raise FormatError, f
        try:
            # check if the desired file exists
            cfile=open(prefix+file_id+'.plot','r')
        except IOError:
            for tmpfile in file_list:
               os.remove(tmpfile)
            msg='Characteristic not computed: '+str(ViewPoint)
            raise ValueError, msg                
        else:
            self.__plot= stat_tool._PlotManager(file_list, prefix+file_id, 
                                                 nb_windows)

##        import os
##        prefix="ftmp"
##        file_created=False
##        file_list=[]
##        while not file_created:
##            try:
##                cfile=open(prefix+'.plot','r')
##            except IOError:
##                # file does not exist
##                for o in range(self.__chmt.NbInt()):
##                    file_list+= [prefix+str(o+1)+extension \
##                                for extension in \
##                                [".plot", ".dat", ".print"]]
##                    file_created=True
##            else:
##                import random
##                prefix+=str(random.randint(1,9))
##        try:
##            self.__chmt.Plot(prefix, Title)
##        except RuntimeError, f:
##            for tmpfile in file_list:
##                os.remove(tmpfile)
##            raise FormatError, f
##        else:                      
##            self.__plot=stat_tool._PlotManager(file_list, 
##                                                prefix+str(variable+1), 1)

    def Save(self, file_name, format="ASCII", overwrite=False):
        """Save HiddenMarkovTree object into a file.
        
        Argument file_name is a string designing the file name and path.
        String argument format must be "ASCII" or "SpreadSheet".
        Boolean argument overwrite is false is the file should not 
        be overwritten."""
        if not overwrite:
            try:
                f=file(file_name, 'r')
            except IOError:
                f=file(file_name, 'w+')
            else:
                msg="File "+file_name+" already exist"
                raise IOError, msg
            f.close()
        if not (string.upper(format)=="ASCII" 
                or string.upper(format)=="SPREADSHEET"):
            msg="unknown file format: "+str(format)
            raise ValueError, msg
        elif (string.upper(format)=="ASCII"):
            try:
                self.__chmt.FileAsciiWrite(file_name)
            except RuntimeError, error:
                raise FormatError, error
        else:
            try:
                self.__chmt.SpreadsheetWrite(file_name)
            except RuntimeError, error:
                raise FormatError, error

        
    def Simulate(self, arg1, arg2, arg3=None):
        """Generate a sample of trees from self.
        
        Usage:  Simulate(sample_size, tree_size, nb_children)
                Simulate(sample_size, Trees)
                Simulate(size_histo, nb_children_histo)"""
        try:
            if (type(arg1)==int) and (type(arg2)==int):
                # Simulate(sample_size, tree_size, nb_children)
                if arg3 is None:
                    raise TypeError, "argument 3 is mandatory in " \
                        "Simulate(int, int, int)"
                chmt_data= self.__chmt.Simulate(arg1, arg2, arg3, True)
            elif type(arg1)==int:
                # Simulate(sample_size, Trees)
                if issubclass(arg2.__class__, trees.Trees):
                    chmt_data= self.__chmt.Simulate(arg1, arg2._ctrees(), True)
                else:
                    raise TypeError, "bad type for argument 2: trees.Trees " \
                                     "expected"
            else:
                # Simulate(size_histo, nb_children_histo)
                if issubclass(arg1.__class__, stat_tool.Histogram):
                    if issubclass(arg2.__class__, stat_tool.Histogram):
                        chmt_data= \
                            self.__chmt.Simulate(arg1._chisto(), arg2._chisto(), 
                                                 True, False)
                    else:
                        raise TypeError, "bad type for argument 2:  " \
                                         "stat_tool.Histogram expected"
                else:
                    raise TypeError, "bad type for argument 1:  " \
                                     "stat_tool.Histogram expected"
        except RuntimeError, error:
                raise FormatError, error
        chmt_data=chmt_data.StateTrees()
        return HiddenMarkovTreeData(chmt_data, self, True)

    def _chmt(self):
        return self.__chmt
    
    def _criteria(self):
        # extract the value of each selection criterion
        disp=self.__chmt.Display(False)
        criteria={}
        names=["AIC", "AICc", "BIC", "BICc", "ICL", "ICLc"]
        for name in names:
            f=disp.find("("+name+"):")
            if (f != -1):
                pos=f+len(name)+3
                i=disp.find("\n", pos)
                try:
                    val=float(disp[pos:i])
                except ValueError:
                    pass
                else:
                    if str(val).upper() != "NAN":
                        criteria[name]=val
        return criteria
    
    def _likelihood(self, trees):
        # compute the likelihood of a given set of trees
        try:
            likelihood=self.__chmt.Likelihood(trees._ctrees())
        except RuntimeError, e:
            raise FormatError, e
        return likelihood
    
    def __str__(self):
        classstr=str(self.__class__)
        res=classstr+": "+str(self._chmt())
        return res


class HiddenMarkovTreeData(trees.Trees):
    """A set of trees associated with one hidden markov tree model."""

    def __init__(self, trees_object, markov=None, aliasing=False, 
                 attribute_names=None):
        ## Initialize a HiddenMarkovTreeData object by copy or 
        ## from (a Trees or a CHmt_data object) and a hidden Markov tree.
        ## Aliasing is used to make an alias between trees_object and
        ## self.__ctrees -> useful for the connection between HiddenMarkovTree
        ## and HiddenMarkovTreeData
        ## The attributes names can be provided
        # and the attribute types should have the possibility to be specified
        if issubclass(trees_object.__class__, HiddenMarkovTreeData):
            # trees_object is supposed to be a HiddenMarkovTreeData object...
            trees.Trees(self, trees_object, attribute_names=attribute_names)
            self.__ctrees=chmt.CHmt_data(arg.__chmtd)
        elif issubclass(trees_object.__class__, chmt.CHmt_data):
            # ...or a CHmt_data object...
            trees.Trees.__init__(self, trees_object,
                                 attribute_names=attribute_names)
            if aliasing:
                self.__ctrees=trees_object
            else:
                self.__ctrees=chmt.CHmt_data(trees_object)
        elif issubclass(trees_object.__class__, HiddenMarkovTree):
            # above line is weird: class must be Trees...
            # ... or a Trees object...
            if trees_object is None:
                raise ValueError, "second argument is mandatory"
            elif issubclass(markov.__class__, HiddenMarkovTree):
                # ... and the second argument must be a HiddenMarkovTree
                # in the latter case
                trees.Trees.__init__(self, trees_object, 
                                     attribute_names=attribute_names)
                self.__ctrees=chmt.CHmt_data(self._ctrees(), markov._chmt())                
            else:
                msg="bad type for second argument: "+str(type(trees_object))
                raise TypeError, msg                
        else:
            msg="bad type for first argument: "+str(type(arg))
            raise TypeError, msg
##        if not (attribute_names is None):
##            self.__attributes=list(attribute_names)
                          
    def ExtractHistogram(self, nature, variable=None, state=None):
        """Extract a frequency distribution from self.
        
        Usage:  ExtractHistogram("Size")
                ExtractHistogram("NbChildren")
                ExtractHistogram("Value", variable)
                ExtractHistogram("VariableName")
                ExtractHistogram("NbZones", variable, value)
                ExtractHistogram("Observation", variable, state)"""
        if type(nature)==str:
            if ((string.upper(nature)=="SIZE") or
                (string.upper(nature)=="NBCHILDREN") or
                (string.upper(nature)=="VALUE")):
                chisto=trees.Trees.ExtractHistogram(self, nature, variable)
                # chisto=self.__ctrees.ExtractValueHistogram(variable)            
            elif ((nature.upper()=="FIRSTOCCURRENCEROOT") or
                  (nature.upper()=="FIRSTOCCURRENCELEAVES") or
                  (nature.upper()=="SOJOURNSIZE") or
                  (nature.upper()=="NBZONES") or
                  (nature.upper()=="NBOCCURRENCES") or
                  (nature.upper()=="OBSERVATION")):
                if variable is None:
                    if ((self.NbVariables()==1) 
                        and (nature.upper()!="OBSERVATION")):
                    # variable argument is not mandatory if there is only
                    # one variable for other features than observation
                        variable=0
                    else:
                        msg='argument 2 is mandatory in ExtractHistogram' + \
                            '(FeatureName, variable, state/value)'
                        raise TypeError, msg
                if state is None:
                    # value argument is mandatory 
                    raise TypeError, 'argument 3 is mandatory in ' \
                        'ExtractHistogram(FeatureName, variable, state/value)'
                elif type(state)!=int:
                    raise TypeError, "bad type for argument 3: " \
                                      "type 'int' expected"
                if nature.upper()=="FIRSTOCCURRENCEROOT":
                    chartype=CharacteristicType.FIRST_OCCURRENCE_ROOT
                elif nature.upper()=="FIRSTOCCURRENCELEAVES":
                    chartype=CharacteristicType.FIRST_OCCURRENCE_LEAVES
                elif nature.upper()=="SOJOURNSIZE":
                    chartype=CharacteristicType.SOJOURN_SIZE
                elif nature.upper()=="NBZONES":
                    chartype=CharacteristicType.NB_ZONES
                elif nature.upper()=="NBOCCURRENCES":
                    chartype=CharacteristicType.NB_OCCURRENCES
                elif nature.upper()=="OBSERVATION":
                    chartype=CharacteristicType.OBSERVATION
                chartype+=0
                try:
                    chisto=self.__ctrees.ExtractValueHistogram(
                        chartype, self._valid_cvariable(variable)+1, state)
                except RuntimeError, error:
                    raise FormatError, error
            else:
                chisto=trees.Trees.ExtractHistogram(self, nature)
            return stat_tool.Histogram(chisto)
        else:
            raise TypeError, "bad type for argument 1: type 'str' expected"

    def Plot(self, ViewPoint="Data", Length=None, BottomDiameter=None, 
             Color=None, DressingFile=None, Title="", variable=0):
        """Graphical output using the Geom 3D viewer for Trees 
           or Gnuplot.py for features.
        
        Usage:  Plot(ViewPoint="Data")
                Plot(ViewPoint="FirstOccurrence", variable=0)
                Plot(ViewPoint="Observation", variable=1)
        Other possible values for ViewPoint: 
            FirstOccurrenceLeaves
            SojournSize
            Counting"""
        if type(ViewPoint)!=str:
            msg='bad type for argument "ViewPoint": '\
              +"type 'str' expected"
            raise TypeError, msg
        if (ViewPoint.upper()!="OBSERVATION"):
            # use Plot for Trees
            trees.Trees.Plot(self, ViewPoint, Length, BottomDiameter,
                             Color, DressingFile, Title, variable)
        else:
            # Plot(ViewPoint="Observation", variable)
            # Graphical output of the features using Gnuplot.py
            import os
            ftype=CharacteristicType.OBSERVATION
            ref_file_id=str(self._valid_cvariable(variable))
            # part of the filename which identifies the graph to be displayed
            prefix="ftmp"
            file_created=False
            file_list=[]
            # find a non existing file name
            while not file_created:
                try:
                    cfile=open(prefix+'1.plot','r')
                except IOError:
                    file_created=True
                else:
                    import random
                    prefix+=str(random.randint(1,9))
            try:
                self.__ctrees.Plot(os.getcwd()+os.sep+prefix, Title)
                # build the list of the created files: 
                for var in [str(self._valid_cvariable(v)+1) 
                            for v in range(self.NbInt())]+["0"]:
                    for char in [str(c) for c in range(5)]+[""]:
                        filename=prefix+var+char
                        try:
                            tmpfile=open(filename+'.plot', 'r')
                        except IOError:
                            pass
                        else:
                            tmpfile.close()
                            # add the .plot and .print files
                            file_list+=[filename+extension
                                for extension in [".plot", ".print"]]
                for var in [str(self._valid_cvariable(v)+1) 
                            for v in range(self.NbInt())]+["0"]:
                    for char in [str(c) for c in range(5)]+[""]:
                        filename=prefix+char+var
                        try:
                            tmpfile=open(filename+'.dat', 'r')
                        except IOError:
                            pass
                        else:
                            # add the .dat file
                            tmpfile.close()
                            file_list+=[filename+".dat"]
            except RuntimeError, f:
                for tmpfile in file_list:
                   os.remove(tmpfile)
                raise FormatError, f
            if (prefix+ref_file_id+'0.plot') in file_list:
                file_id=ref_file_id+"0"
            else:
                file_id=ref_file_id
            try:
                # check if the desired file exists
                cfile=open(prefix+file_id+'.plot','r')
            except IOError:
                for tmpfile in file_list:
                   os.remove(tmpfile)
                if variable==0:
                    msg='The hidden variable is no output process'
                else:
                    msg='Observation distribution not available: '+str(ViewPoint)
                raise ValueError, msg                
            else:
                # plot the files using the PlotManager class of stat_tool
                nb_windows= \
                   self.__ctrees.NbValues(0)
                self.__plot=stat_tool._PlotManager(file_list, prefix+file_id, 
                                                    nb_windows)
                                                     
    def _state_marginal_distribution(self):
        # compute the (empirical) marginal distribution of the hidden states
        s=ctrees.CTrees.Display(self._ctrees(), False)
        i=s.find("state histogram - sample size", 0)
        msg="Could not find the (empirical) marginal distribution" \
            " of the hidden states"
        if i == -1:
            raise Warning, msg
        i=s.find("state histogram", i+1)
        if i == -1:
            raise Warning, msg
        nb_states=self.__ctrees.NbValues(0)
        res=range(nb_states)
        total=0.
        for v in range(nb_states):
            i=s.find(str(v), i+1)
            e=s.find("\n", i+1)
            try:
                res[v]=int(s[i+len(str(v)):e])
            except ValueError:
                raise Warning, msg
            total+=res[v]
            i=e
        for v in range(nb_states):
            res[v]=res[v] / total
        return res
    
    def __str__(self):
        classstr=str(self.__class__)
        res=classstr+": "+str(self._ctrees())
        return res


