# -*- coding: utf-8 -*-
"""Hidden Markov tree models"""
import string
import openalea.stat_tool.error as check_error
import openalea.tree_statistic._errors as _errors
import openalea.stat_tool as stat_tool, openalea.tree_statistic.trees as trees
import openalea.stat_tool.error as check_error
import _hmt

VariableType=stat_tool.VariableTypeBis
StatTreeError = _errors.StatTreeError
CharacteristicType=trees.CharacteristicType
EntropyAlgorithm=_hmt.EntropyAlgorithm
VariableTypeDict=VariableType.values

from openalea.stat_tool import interface
interface.extend_class(_hmt.CiHmot, interface.StatInterface)
interface.extend_class(_hmt.CHmt_data, interface.StatInterface)

class HiddenMarkovIndOutTree:
    """An implementation of the hidden Markov out-trees with conditionally 
    independent children states given their parent."""

    def __init__(self, arg, aliasing=False):
        """Initialize a Hmt by copy or by reading into a file.

            Usage:  H=HiddenMarkovIndOutTree("file_name")
                    H=HiddenMarkovIndOutTree(HiddenMarkovIndOutTree)."""
        ## Aliasing is used to make an alias between argument and
        ## self.__chmt -> useful for the connection between HiddenMarkovIndOutTree
        ## and HiddenMarkovTreeData
        if type(arg)==str:
            # arg is expectedly a file name...
            self.__chmt=_hmt.HmtAsciiRead(arg)
            nbvariables=self.__chmt.NbInt()+self.__chmt.NbFloat()
            self._attributes=["Variable " + str(i) for i in range(nbvariables)]
        elif issubclass(arg.__class__, HiddenMarkovIndOutTree):
            # ... or a HiddenMarkovIndOutTree object...
            if aliasing:
                self.__chmt=arg.__chmt
            else:
                self.__chmt=_hmt.CiHmot(arg.__chmt)
            self._attributes=list(arg._attributes)
        elif issubclass(arg.__class__, _hmt.CiHmot):
            # ... or a CiHmot object
            if aliasing:
                self.__chmt=arg
            else:
                self.__chmt=_hmt.CiHmot(arg)
            nbvariables=self.__chmt.NbInt()+self.__chmt.NbFloat()
            self._attributes=["Variable " + str(i) for i in range(nbvariables)]
        else:
            msg="bad argument type: "+str(type(arg))
            raise TypeError, msg

    def Display(self, ViewPoint=None, Detail=None, TreeId=None, 
                NbStateTrees=2, StateTrees="GeneralizedViterbi",
                Entropy="UPWARD", RootVertex=None):
        """Display HiddenMarkovIndOutTree object.
        
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
                if (type(TreeId) != int):
                    msg='bad type for argument "TreeId": '\
                        +str(type(TreeId))
                    raise TypeError, msg
                if (type(NbStateTrees)!=int):
                    msg='bad type for argument "NbStateTrees": '\
                        +str(type(NbStateTrees))
                    raise TypeError, msg
                if not(RootVertex is None) and (type(RootVertex)!=int):
                    msg='bad type for argument "RootVertex": '\
                        +str(type(RootVertex))
                    raise TypeError, msg
                elif RootVertex is None:
                    rv = stat_tool.I_DEFAULT
                else:
                    rv = RootVertex
                if (type(Entropy)!=str):
                    msg='bad type for argument "Entropy": '\
                        +str(type(Entropy))
                    raise TypeError, msg
                elif (Entropy.upper()=="UPWARD"):
                    entropy_algo=EntropyAlgorithm.UPWARD
                elif (Entropy.upper()=="DOWNWARD"):
                    entropy_algo=EntropyAlgorithm.DOWNWARD
                else:
                    msg="Bad value for 'Entropy' argument:"+str(Entropy) \
                    +" - expecting 'UPWARD' or 'DOWNWARD'"
                    raise ValueError, msg    
                res = self.__chmt.StateProfile(0, NbStateTrees,
                                                TreeId, entropy_algo,
                                                rv)
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

    def EntropyComputation(self, TreeId=None, UpwardEntropy=False):
        """
        Compute state entropy, sum of upward entropies and sum of
        marginal entropies for a given tree or the whole set of trees.
        """
        check_error.CheckType([UpwardEntropy], [bool])
        if (UpwardEntropy):
            if TreeId is None:
                res = self.__chmt.UpwardEntropyComputation()
            else:
                check_error.CheckType([TreeId], [int])
                res = self.__chmt.UpwardEntropyComputation(TreeId)
        else:
            if TreeId is None:
                res = self.__chmt.EntropyComputation()
            else:
                check_error.CheckType([TreeId], [int])
                res = self.__chmt.EntropyComputation(TreeId)
        return res

    def Extract(self, nature, variable, value):
        """Extract a distribution from the HiddenMarkovIndOutTree.

        :Parameters:

        * `nature` (str) - 'FirstOccurrenceRoot' or 'FirstOccurrenceLeaves' or
                           'SojournSize' or 'nbzones' or 'NbOccurrences'
                           'Observation'
        * `variable` (int) - variable index,
        * `value` (int) - value of variable,

        :Returns:

            If the characteristic is available, its distribution for
            the given variable and value of this variable is returned.

        :Examples:

        .. doctest::
            :options: +SKIP

            >>> 
            >>> self.Extract('FirstOccurrenceRoot', 0, 1)
            >>> self.Extract('Observation', 1, 1)

        .. seealso::

        :func:`~openalea.tree_statistic.hmt.HiddenMarkovTreeData.ExtractHistogram
        :func:`~openalea.tree_statistic.hmt.CharacteristicType
        """

        if type(nature)==str:
            if type(variable)!=int:
                raise TypeError, "bad type for argument variable: " \
                                    "type 'int' expected"
            if type(value)!=int:
                raise TypeError, "bad type for argument value: " \
                                    "type 'int' expected"
            if nature.upper()=="FIRSTOCCURRENCEROOT":
                chartype = CharacteristicType.FIRST_OCCURRENCE_ROOT
            elif nature.upper()=="FIRSTOCCURRENCELEAVES":
                chartype = CharacteristicType.FIRST_OCCURRENCE_LEAVES
            elif nature.upper()=="SOJOURNSIZE":
                chartype = CharacteristicType.SOJOURN_SIZE
            elif nature.upper()=="NBZONES":
                chartype = CharacteristicType.NB_ZONES
            elif nature.upper()=="NBOCCURRENCES":
                chartype = CharacteristicType.NB_OCCURRENCES
            elif nature.upper()=="OBSERVATION":
                chartype = CharacteristicType.OBSERVATION
                chartype+=0
            else:
                msg='bad value for argument "Nature": '+ nature
                raise ValueError, msg
            distrib = self.__chmt.Extract(chartype, variable+1, value)
        else:
            raise TypeError, "bad type for argument 1: type 'str' expected"
        return distrib

    def ExtractData(self):
        """Extract the 'data' part of the HiddenMarkovIndOutTree."""
        chmt_data = self.__chmt.ExtractData()
        chmt_data = chmt_data.StateTrees()
        return HiddenMarkovTreeData(chmt_data, markov=None, aliasing=True)

    def ExtractPlot(self, TreeId, ViewPoint, Entropy="UPWARD"):
        """Extract a tree with state or entropy profile.
        
        Usage:  ExtractPlot(TreeId, ViewPoint="StateProfile")"""
        if (TreeId is None):
            raise TypeError, "Argument 'TreeId' is mandatory in "+\
                    'ExtractPlot(TreeId, ViewPoint="StateProfile")'
        if (type(Entropy)!=str):
            msg='bad type for argument "Entropy": '\
                +str(type(Entropy))
            raise TypeError, msg
        elif (Entropy.upper()=="UPWARD"):
            entropy_algo=EntropyAlgorithm.UPWARD
        elif (Entropy.upper()=="DOWNWARD"):
            entropy_algo=EntropyAlgorithm.DOWNWARD                    
        res = self.__chmt.StateProfile(0, 2, TreeId, entropy_algo)
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

    def NbStates(self):
        """Return the number of hidden states of the HiddenMarkovIndOutTree."""
        return self.__chmt.NbStates()
    
    def Plot(self, ViewPoint="Observation", variable=0, Title="", 
             TreeId=0, EndVertex=None, Entropy="UPWARD"):
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
            ftype = CharacteristicType.FIRST_OCCURRENCE_ROOT
        elif ViewPoint.upper()=="FIRSTOCCURRENCELEAVES":
            ftype = CharacteristicType.FIRST_OCCURRENCE_LEAVES
        elif ViewPoint.upper()=="SOJOURNSIZE":
            ftype = CharacteristicType.SOJOURN_SIZE
        elif ViewPoint.upper()=="COUNTING":
            ftype = CharacteristicType.NB_ZONES
        elif ViewPoint.upper()=="OBSERVATION":
            ftype=0 #CharacteristicType.OBSERVATION
            is_observation=True
        elif ViewPoint.upper()!="STATEPROFILE":
            msg='bad value for argument "ViewPoint": '+ViewPoint
            raise ValueError, msg
        if ViewPoint.upper()!="STATEPROFILE":
            try:
                w = self.__chmt.NbValues(variable)
            except StatTreeError:
                msg = "variable index out of range: "+str(variable)
                raise StatTreeError(msg)
            ref_file_id=str(variable+1)
            # part of the filename which identifies the graph to be displayed
            file_id=ref_file_id
            # in the special case of observation distributions,
            # the suffix can be "0" or void, depending on
            # the number of values
            if not(is_observation):
                file_id+=str(ftype+2)
            elif not(self.__chmt.IsParametric(variable)):
                file_id+=str(0)
            self.__chmt.plot(Title=Title, Suffix=file_id)
        else: # ViewPoint.upper()=="STATEPROFILE"
            if (type(TreeId)!=int):
                raise TypeError, "bad tree index type: "+str(type(TreeId))
            elif (TreeId < 0):
                raise IndexError, "bad tree index value: "+str(variable)
            elif (type(EndVertex)!=int):
                msg="bad type for vertex identifier: " + \
                         str(type(EndVertex))
                raise TypeError, msg
            if (type(Entropy)!=str):
                msg='bad type for argument "Entropy": '\
                    +str(type(Entropy))
                raise TypeError, msg
            elif (Entropy.upper()=="UPWARD"):
                entropy_algo=EntropyAlgorithm.UPWARD
            elif (Entropy.upper()=="DOWNWARD"):
                entropy_algo=EntropyAlgorithm.DOWNWARD
            else:
                msg="Bad value for 'Entropy' argument:"+str(Entropy) \
                +" - expecting 'UPWARD' or 'DOWNWARD'"
                raise ValueError, msg    
            ref_file_id=""
            file_id=ref_file_id
            self.__chmt.plot(Title=Title, Suffix=file_id,
                                ViewPoint="StateProfile",
                                Params=(TreeId, EndVertex, entropy_algo))

    def Print(self, ViewPoint="Observation", variable=0, Title=""):
        """Graphical output into a postscript file using Gnuplot.py.        
        
        Usage:  Print(ViewPoint="FirstOccurrenceRoot", variable=0)
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
            ftype = CharacteristicType.FIRST_OCCURRENCE_ROOT
        elif ViewPoint.upper()=="FIRSTOCCURRENCELEAVES":
            ftype = CharacteristicType.FIRST_OCCURRENCE_LEAVES
        elif ViewPoint.upper()=="SOJOURNSIZE":
            ftype = CharacteristicType.SOJOURN_SIZE
        elif ViewPoint.upper()=="COUNTING":
            ftype = CharacteristicType.NB_ZONES
        elif ViewPoint.upper()=="OBSERVATION":
            ftype=0 #CharacteristicType.OBSERVATION
            is_observation=True
        else:
            msg='bad value for argument "ViewPoint": '+ViewPoint
            raise ValueError, msg
        try:
            w = self.__chmt.NbValues(variable)
        except StatTreeError:
            msg = "variable index out of range: "+str(variable)
            raise IndexError, msg          
        # part of the filename which identifies the involved variable
        ref_file_id = str(variable+1)
        # part of the filename which identifies the graph to be displayed
        file_id=ref_file_id
        if not(is_observation):
            file_id+=str(ftype+2)
        elif not(self.__chmt.IsParametric(variable)):
            file_id+=str(0)
        # in the special case of observation distributions,
        # the suffix can be "0" or void, depending on the number of values
        self.__chmt.plot_print(Title=Title, Suffix=file_id)


    def Save(self, file_name, format="ASCII", overwrite=False):
        """Save HiddenMarkovIndOutTree object into a file.
        
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
                msg="File "+file_name+" already exists"
                raise IOError, msg
            f.close()
        import string
        if not (string.upper(format)=="ASCII" 
                or string.upper(format)=="SPREADSHEET"):
            msg="unknown file format: "+str(format)
            raise ValueError, msg
        elif (string.upper(format)=="ASCII"):
            self.__chmt.FileAsciiWrite(file_name)
        else:
            self.__chmt.SpreadsheetWrite(file_name)


    def Simulate(self, arg1, arg2=None, arg3=None):
        """Generate a sample of trees from self.

        Usage:  Simulate(sample_size, tree_size, nb_children)
                Simulate(sample_size, Trees)
                Simulate(size_histo, nb_children_histo)
                Simulate(Trees)"""
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
        elif issubclass(arg1.__class__, stat_tool.histogram._DiscreteDistributionData):
            # Simulate(size_histo, nb_children_histo)
            expected_class = stat_tool.histogram._DiscreteDistributionData
            if issubclass(arg1.__class__, expected_class):
                if issubclass(arg2.__class__, expected_class):
                    chmt_data= \
                        self.__chmt.Simulate(arg1, arg2, True, False)
                else:
                    msg = "bad type for argument 2:  " +  \
                          str(expected_class) + "expected"
                    raise TypeError, msg
            else:
                msg = "bad type for argument 1:  " +  \
                        str(expected_class) + "expected"
                raise TypeError, msg
        else:
            # Simulate(Trees)
            expected_class = trees.Trees
            msg = "bad type for argument 1:  " +  \
                       str(expected_class) + "expected"
            if issubclass(arg1.__class__, expected_class):
                chmt_data = self.__chmt.Simulate(arg1._ctrees(), True)
            else:
                raise TypeError, msg
        chmt_data = chmt_data.StateTrees()
        return HiddenMarkovTreeData(chmt_data, self, True)

    def StatePermutation(self, perm):
        """Permutation of the states of self.
        perm[i]==j means that current state i will become new state j.

        Usage:  StatePermutation(list)"""
        self.__chmt.StatePermutation(perm)

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
        likelihood = self.__chmt.Likelihood(trees._ctrees())
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
        ## self.__ctrees -> useful for the connection between HiddenMarkovIndOutTree
        ## and HiddenMarkovTreeData
        ## The attributes names can be provided
        # and the attribute types should have the possibility to be specified
        if issubclass(trees_object.__class__, HiddenMarkovTreeData):
            # trees_object is supposed to be a HiddenMarkovTreeData object...
            trees.Trees.__init__(self, trees_object, 
                                 attribute_names=attribute_names)
            self.__ctrees=_hmt.CHmt_data(trees_object._ctrees(), True)
        elif issubclass(trees_object.__class__, _hmt.CHmt_data):
            # ...or a CHmt_data object...
            trees.Trees.__init__(self, trees_object,
                                 attribute_names=attribute_names)
            if aliasing:
                self.__ctrees=trees_object
            else:
                self.__ctrees=_hmt.CHmt_data(trees_object, True)
        elif issubclass(trees_object.__class__, trees.Trees):
            # ... or a Trees object...
            if trees_object is None:
                raise ValueError, "second argument is mandatory"
            elif issubclass(markov.__class__, HiddenMarkovIndOutTree):
                # ... and the second argument must be a HiddenMarkovIndOutTree
                # in the latter case
                trees.Trees.__init__(self, trees_object, 
                                     attribute_names=attribute_names)
                self.__ctrees=_hmt.CHmt_data(super(HiddenMarkovTreeData, self)._ctrees(), markov._chmt())
            else:
                msg="bad type for second argument: "+str(type(trees_object))
                raise TypeError, msg                
        else:
            msg="bad type for first argument: "+str(type(trees_object))
            raise TypeError, msg
##        if not (attribute_names is None):
##            self.__attributes=list(attribute_names)

    def ExtractMarkov(self):
        """Extract the 'model' part of the HiddenMarkovTreeData."""
        chmt = self.__ctrees.ExtractMarkov()
        return HiddenMarkovIndOutTree(chmt, aliasing=True)

    def ExtractHistogram(self, nature, variable=None, state=None):
        """Extract a frequency distribution from self.
        
        Usage:  ExtractHistogram("Size")
                ExtractHistogram("NbChildren")
                ExtractHistogram("Value", variable)
                ExtractHistogram("VariableName")
                ExtractHistogram("NbZones", variable, value)
                ExtractHistogram("Observation", variable, state)"""
        check_error.CheckType([nature], [str])
        if ((string.upper(nature)=="SIZE") or
            (string.upper(nature)=="NBCHILDREN")):
            chisto = trees.Trees.ExtractHistogram(self, nature, variable)
            # chisto=self.__ctrees.ExtractValueHistogram(variable)
        elif (string.upper(nature)=="VALUE"):
            # Extract marginal histogram with mixture of observation distributions
            check_error.CheckType([variable], [int])
            if (variable < 1) or (variable > self.NbVariables()):
                msg = "bad variable: " + str(variable)
                raise IndexError, msg
            chisto = self.__ctrees.ExtractMarginal(variable+1)
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
                    variable = 0
                else:
                    msg = 'argument 2 is mandatory in ExtractHistogram' + \
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
                chartype = CharacteristicType.FIRST_OCCURRENCE_ROOT
            elif nature.upper()=="FIRSTOCCURRENCELEAVES":
                chartype = CharacteristicType.FIRST_OCCURRENCE_LEAVES
            elif nature.upper()=="SOJOURNSIZE":
                chartype = CharacteristicType.SOJOURN_SIZE
            elif nature.upper()=="NBZONES":
                chartype = CharacteristicType.NB_ZONES
            elif nature.upper()=="NBOCCURRENCES":
                chartype = CharacteristicType.NB_OCCURRENCES
            elif nature.upper()=="OBSERVATION":
                chartype = CharacteristicType.OBSERVATION
            chartype+=0
            chisto = self.__ctrees.ExtractValueHistogram(
                        chartype, self._valid_cvariable(variable)+1, state)
        else:
            chisto = trees.Trees.ExtractHistogram(self, nature)
        return chisto

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
            ftype = CharacteristicType.OBSERVATION
            if ((self.Type(variable) == VariableType.INT_VALUE) or
                (self.Type(variable) == VariableType.REAL_VALUE)):
                cvariable = self._valid_cvariable(variable)
            else:
                msg = "Bad variable type for variable ", variable, ": " + \
                      str(self.Type(variable))
                raise TypeError, msg
            ref_file_id = str(cvariable)
            # if the number of values is low for variable, 
            # each characteristic is drawn into a file with "0" suffix
            if (self.__ctrees.IsCharacteristic(cvariable,
                                               CharacteristicType.NB_ZONES)):
                file_id = ref_file_id + "0"
            else:
                file_id = ref_file_id
            # try:
            #     param = False
            param = self.__ctrees.IsParametric(cvariable)
            # except RuntimeError:
            #     pass
            # else:
            if param:
                file_id = ref_file_id
            # part of the filename that identifies the graph to be displayed
            self.__ctrees.plot(Title=Title, Suffix=file_id)

    def _ctrees(self):
        return self.__ctrees

    def _state_marginal_distribution(self):
        # compute the (empirical) marginal distribution of the hidden states
        s=trees.Trees._ctrees_display(self)
        msg="Could not find the (empirical) marginal distribution" \
            " of the hidden states"
        i = s.find("VARIABLE 1 : STATE", 0)
        if i == -1:
            raise Warning, msg
        i = s.find("state marginal frequency distribution - sample size", i+1)
        if i == -1:
            raise Warning, msg
        i=s.find("| frequency", i+1)
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


