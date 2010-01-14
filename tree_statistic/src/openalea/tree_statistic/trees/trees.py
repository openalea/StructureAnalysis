# -*- coding: utf-8 -*-
"""Trees with integral and floating attributes and Tree-related objects"""
__revision__ = "$Id$"

import string
import openalea.stat_tool as stat_tool
import openalea.tree_statistic.int_fl_containers as int_fl_containers
import ctree, ctrees

from openalea.stat_tool import interface
interface.extend_class(ctrees.CTrees, interface.StatInterface)

# _PlotManager=stat_tool.stat_tool._PlotManager
I_DEFAULT_TREE_SIZE = ctree.I_DEFAULT_TREE_SIZE()
I_DEFAULT_TREE_SIZE = ctree.I_DEFAULT_TREE_SIZE()
I_DEFAULT_TREE_DEPTH = ctree.I_DEFAULT_TREE_DEPTH()
VariableType = stat_tool.VariableTypeBis
FormatError = stat_tool.FormatError
CharacteristicType = ctree.Characteristic

VariableTypeDict = VariableType.values

class TreeValue:
    """A class for handling the attributes of class Tree, like a list of int
    and double with type checking."""

    def __init__(self, list_of_values):
        """
        Initialize a container of integral and floating values.

        :Parameters:
          `list_of_values` (list) - determines the type and the actual value \
            of each variable.

        :Examples:
        
        .. doctest::
        
            >>> v = TreeValue([0, 0., 1])
        """
        self.__values = []
        self.__types = []
        self.__nb_integral = 0
        self.__nb_float = 0
        nb_variables = 0
        if ((not hasattr(list_of_values, "__getitem__")) and 
            (not issubclass(list_of_values.__class__, TreeValue))):
            msg="bad type for attribute list: "+str(type(list_of_values))
            raise TypeError, msg
        for index in range(len(list_of_values)):
            o = list_of_values[index]
            if issubclass(o.__class__, VariableType):
                self.__types.append(o)
                if (o == VariableType.REAL_VALUE):
                    val = 0.
                else:
                    val = 0
                self.__values.append(val)
                nb_variables += 1
            elif type(o) == int:
                self.__types.append(VariableType.INT_VALUE)
                self.__values.append(o)
                self.__nb_integral += 1
                nb_variables += 1
            elif type(o) == float:
                self.__types.append(VariableType.REAL_VALUE)
                self.__values.append(o)
                self.__nb_float += 1
                nb_variables += 1
            else:
                s = "element %d of the list of values must be of type int or "\
                "double"%nb_variables
                raise TypeError, s
                # STAT_error[STATR_VARIABLE_TYPE];

    def NbInt(self):
        """Return the number of integral values."""
        return self.__nb_integral

    def NbFloat(self):
        """Return the number of floating values."""
        return self.__nb_float

    def Types(self):
        """Return the list of the container types."""
        return self.__types

    def Type(self, index):
        """
        Return the type of one given variable.

        :Parameters:
          `index` (int) - refers to the concerned variable. \
            `index` must be in the range [0, len(self)].

        :Returns:
            If `index` must be in the range [0, len(self)], a variable type
            among :ref:`openalea.stat_tool._stat_tool.VariableTypeBis` is returned.
            
        :Examples:

        .. doctest::
        
            >>> v = TreeValue([0, 0., 1])
            >>> v.Type(0) == VariableTypeBis.INT_VALUE
        """
        return self.__types[index]

    def Values(self):
        """Return the values of the container"""
        return self.__values

    def __getitem__(self, index):
        return self.__values[index]

    def __setitem__(self, index, valeur):
        if type(valeur) == int:
            if self.__types[index] == VariableType.INT_VALUE:
                # or other integral types
                self.__values[index] = valeur
            elif self.__types[index] == VariableType.REAL_VALUE:
                self.__values[index] = valeur+0.
                # conversion from int to double is allowed
            else:
                raise TypeError, "expected type: INT_VALUE or REAL_VALUE"
                # should not happen if self has been created properly
        elif type(valeur) == float:
            if self.__types[index] == VariableType.REAL_VALUE:
                self.__values[index] = valeur
            else:
                raise TypeError, "expected type: INT_VALUE"
        else:
            raise TypeError, "expected type: INT_VALUE or REAL_VALUE"

    def __len__(self):
        return len(self.__values)

    def __str__(self):
        return str(self.__values)



class Tree:
    """An implementation of trees with integral and floating attributes."""

    def __init__(self, arg, arg2=None, 
                 attribute_names=None, attribute_def=None, scale=None):
        """Define a Tree from a list of values and a tree object.

        Initialize Tree object from either:
        - a list of integral and floating values or a TreeValue object
          and a TreeStructure object;
        - a Tree object;
        - a MTG file, a filter on the vertices, a list of attribute names,
          a list of attribute functions and the considered scale.
        Example: T = Tree([0, 1.], TreeStructure(TS))"""
        self.__mtg_to_tree_vid = None
        self.__tree_to_mtg_vid = None
        self.__mtg_tid = None
        if attribute_names is None:
            attributes=[]
        else:
            attributes = list(attribute_names)
        if issubclass(arg.__class__, Tree):
            # arg is supposed to be a Tree...
            self.__ctree = ctree.CTree(arg.__ctree)
            self.__types = list(arg.__types)
            if len(attributes) > 0:
                attributes = list(attribute_names)
            else:
                attributes = list(arg.__attributes)
            self._copy_vid_conversion(arg.__mtg_to_tree_vid, 
                                      arg.__tree_to_mtg_vid)
            self.__mtg_tid = arg.__mtg_tid
        elif type(arg)==str:
            # ... or the name of a MTG...
            trees_object = Trees(arg, arg2, attribute_names, attribute_def, scale)
            tree_object = trees_object.Tree(0)
            self.__ctree = ctree.CTree(tree_object.__ctree)
            self.__types = list(tree_object.__types)
            self._copy_vid_conversion(trees_object.TreeVertexId(0),
                                      trees_object.MTGVertexId(0))
            #self.__mtg_tid = trees_object.MTGComponentRoot(0)
            if len(attributes) == 0:
                self.__attributes = list(tree_object.Attributes())
                attributes = list(self.__attributes)
            if trees_object.NbTrees() > 1:
                msg="MTG "+str(arg)+" corresponds to a forest, not to a single"\
                    + " tree structure.\n Use constructor trees.Trees()" \
                    + " instead of trees.Tree()"
                import sys
                sys.stderr.write(msg)
                raise Tree, self
        else:
            if not hasattr(arg, "__getitem__"):
                raise TypeError, "argument 1 must have a __getitem__ method"
            if issubclass(arg[0].__class__, VariableType):
                #... or a list of types
                self.__types = list(arg)
                default_value = TreeValue(arg)
            else:
                #... or a list of values
                default_value = TreeValue(arg)
                self.__types = list(default_value.Types())
            self.__attributes=[]
            # arg2 is supposed to be a tree structure...
            i = int_fl_containers.Int_fl_container(default_value)
            if arg2 is None:
                arg2 = TreeStructure()
            elif issubclass(arg2.__class__, ctree.CTree):
                # ... or a ctree.CTree object
                self.__ctree = ctree.CTree(arg2)
            else:
                self.__ctree = ctree.CTree(arg2._tree(), i)
        if len(attributes) == 0:
            # Default attribute name : "Variable0", etc.
            for var in range(self.NbVariables()):
                attributes.append("Variable"+str(var))
        self.__attributes = list(attributes)
        
    def Attributes(self):
        """
        Return the name of the tree attributes.

        :Returns:
            A list of `str` is returned
        """
       
        return list(self.__attributes)

    def Depth(self):
        """Return the tree depth."""
        return self.__ctree.Depth()

    def Display(self, vids=True, attributes=True, mtg_vids=False):
        """
        Display the tree with or without the vertex identifiers (vids)
        and the attributes.

        :Parameters:
          * `vids` (bool) - If attributes is True, the vids are displayed iif
            vids is True,
          * `attributes` (bool) - the attributes are displayed (after the vids) \
            iif `attributes` is True. If attributes is False, only the vids \
            are displayed,
          * `mtg_vids` (bool) -  if mtg_vids and vids are both True, the MTG \
            vids are displayed instead of the tree vids.
        """
        print self._display(self.Root(), vids, attributes, mtg_vids)

    def EdgeType(self, parent, child):
        """
        Return the type of one given edge (parent, child).

        :Parameters:
          `parent` (int) and `child` (int) - two vertex identifiers \
            that define an edge of self.

        :Returns:
           If (parent, child) defines a valid edge of self, \
           its type is returned (str '+' or '<').
        """
        self.__valid_edge(parent, child)
        btype = self.__ctree.EdgeType(parent, child)
        if btype:
            return "<"
        else:
            return "+"

    def Get(self, vid):
        """
        Return the values of the attributes at a given vertex.

        :Parameters:
          `vid` (int) - a valid vertex identifier (vid).

        :Returns:
            If `vid` defines a valid vertex of self, a :ref:`openalea.tree_statistic.trees.TreeValue` \
            object is returned.
            
        :Examples:

        .. doctest::
            :options: +SKIP

            >>> T.Get(T.Root())
        """
        self.__valid_vid(vid)
        i = self.__ctree.Get(vid)
        values = []
        nb_integral = 0
        nb_float = 0
        for v in range(i.NbInt()+i.NbFloat()):
            if self.__types[v] == VariableType.REAL_VALUE:
                values.append(i.Double(nb_float))
                nb_float += 1
            else:
                values.append(i.Int(nb_integral))
                nb_integral += 1
        res = TreeValue(values)
        return res

    def MTGVertex(self, treevid=None):
        """
        Return the MTG vertex identifier (vid) of a Tree vertex

        :Parameters:
          `treevid` (int) - a valid vertex identifier (vid).

        :Returns:
            If `treevid` defines a valid vertex of self, and that self
            was obtained from a MTG, a MTG vid (int) is returned.
        """
        if self.__tree_to_mtg_vid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if treevid is None:
            return dict(self.__tree_to_mtg_vid)
        else:
            return self.__tree_to_mtg_vid[treevid]
        Argument 

    def NbFloat(self):
        """Return the number of variables with floating type."""
        return self.__ctree.NbFloat()

    def NbInt(self):
        """Return the number of variables with integer type."""
        return self.__ctree.NbInt()

    def NbVariables(self):
        """Return the number of variables (attributes) of a tree."""
        return len(self.__types)

    def Put(self, vid, value):
        """
        Set the attribute values of a given vertex.

        :Parameters:
          * `vid` (int) - a valid vertex identifier (vid),
          * `value` (:ref:`TreeValue` object or list of int and float) - \
            future value of the attributes of self at vertex `vid`

        :Examples:

        .. doctest::
            :options: +SKIP

            >>> T.Put(R.Root(), [0])
        """
        self.__valid_vid(vid)
        tvalue = self.__valid_value(value)
        self.__ctree.Put(vid, tvalue)

    def Round(self, ndigits=0):
        """Round each floating variable to a given precision in decimal digits"""
        res=Tree(self)
        vertex_it= res.__ctree.Vertices()
        types=res.Types()
        ftypes=[]
        for t in range(len(types)):
            if types[t]==VariableType.REAL_VALUE:
                ftypes+=[t]
        try:
            while 1:
                current_vertex=vertex_it.next()
                val=res.Get(current_vertex)
                for v in ftypes:
                    val[v]=round(val[v], ndigits)
                val=res.Put(current_vertex, val)
        except StopIteration: pass
        return res
    
    def Root(self):
        """Return the root vertex identifier (vid)."""
        return self.__ctree.Root()

    def Save(self, file_name, overwrite=False, variable_names=None):
        """Save tree into a file as a MTG."""
        if variable_names is None:
            variable_names=self.__attributes
        elif len(variable_names)!=self.NbVariables():
            msg="bad number of variable names: "+str(len(variable_names)) \
                +"; "+str(self.NbVariables())+" name(s) expected"
            raise IndexError, msg
        for var in range(len(variable_names)):
            if type(variable_names[var])!=str:
                msg="bad type for name of variable "+str(var)
                raise TypeError, msg
            elif variable_names[var].find(" ")!=-1:
                msg="name of variable "+str(var)+" should not contain "\
                     +"any space character"
                raise ValueError, msg
        if not overwrite:
            try:
                f=file(file_name, 'r')
            except IOError:
                pass # f=file(file_name, 'w+')
            else:
                f.close()
                msg="File "+file_name+" already exists"
                raise IOError, msg
                                
        f=file(file_name, 'w+')
        f.write(self.__mtg_header(variable_names))
        f.write(self._mtg_write())
        f.close()
        
    def SelectSubTree(self, vid, keep=True):
        """Select and return a subtree.

        Argument vid, which must be a valid vertex identifier, is the root
        of the selected subtree. Boolean argument keep is True if the
        subtree must be kept, False if it must be pruned from tree."""
        # build the subtree
        self.__valid_vid(vid)
        csubtree=ctree.SelectSubTree(self.__ctree, vid, keep)
        res=Tree(self.Get(self.Root()))
        res.__ctree=ctree.CTree(csubtree)
        # build the attributes
        res.__attributes=list(self.__attributes)
        # build the dictionnaries corresponding to the tree -> MTG
        # and MTG -> tree vid conversions
        if not(self.__mtg_to_tree_vid is None):
            selfvertices=[vid]
            resvertices=[res.Root()]
            res.__mtg_to_tree_vid={}
            res.__tree_to_mtg_vid={}
            def SubTreeTraversal(self, svertices, rvertices,
                                 mtg_to_tree_vid, tree_to_mtg_vid):
                if len(svertices) > 0:
                    # current tree vid
                    v=svertices.pop(0)
                    # corresponding vid for res tree 
                    r=rvertices.pop(0)
                    tree_to_mtg_vid[r]=self.__tree_to_mtg_vid[v]
                    mtg_to_tree_vid[tree_to_mtg_vid[r]]=r
                    svertices+=self.__ctree.Children(v)
                    rvertices+=res.__ctree.Children(r)
                    SubTreeTraversal(self, svertices, rvertices,
                                     mtg_to_tree_vid, tree_to_mtg_vid)
                
                
            SubTreeTraversal(self, selfvertices, resvertices, 
                             res.__mtg_to_tree_vid,
                             res.__tree_to_mtg_vid)
            res.__mtg_tid=self.__mtg_tid
        return res

    def Simulate(self, distributions):
        """Simulate the integer attributes independently.

        Argument distributions is a list of distributions of the same
        length as self.NbInt()
        """
        self.__ctree.IidSimulation(distributions)

    def Size(self):
        """Return the number of vertices."""
        return self.__ctree.Size()

    def TreeVertex(self, mtgvid=None):
        """Return the tree vid of a MTG vertex"""
        if self.__mtg_to_tree_vid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if mtgvid is None:
            return dict(self.__mtg_to_tree_vid)
        else:
            return self.__mtg_to_tree_vid[mtgvid]

    def Type(self, index):
        """Return the type of one given variable.

        Argument index, which must be in the range [0, NbVariables()],
        refers to the concerned variable."""
        return self.__types[index]

    def Types(self):
        """Return the type of each tree attribute."""
        return self.__types

    def _ctree(self):
        return self.__ctree

    def _copy_mtg_tid(self, vid):
        # copy the vid in the MTG corresponding to the root of self
        self.__mtg_tid=vid
        
    def _copy_vid_conversion(self, source_mtg2treedict, source_tree2mtgdict):
        # copy the dictionaries corresponding to the tree -> MTG
        # and MTG -> tree vid conversion from source to self
            if source_mtg2treedict is None:
                self.__mtg_to_tree_vid=None
            else:
                self.__mtg_to_tree_vid=dict(source_mtg2treedict)
            if source_tree2mtgdict is None:
                self.__tree_to_mtg_vid=None
            else:
                self.__tree_to_mtg_vid=dict(source_tree2mtgdict)

    def _display(self, key, vids=True, attributes=True, mtg_vids=False):
        # display the tree, starting from key, and adds a header before
        
        # determine the smoothed probabilities 
        # and quantities to be rounded for display
        smoothed=[]
        rounded=[]
        for index in range(len(self.__attributes)):
            if (self.Type(index)==VariableType.REAL_VALUE):
                if (self.__attributes[index].find("State ")!=-1):
                    smoothed.append(index)
                elif (self.__attributes[index].find("Entropy")!=-1):
                    rounded.append(index)
        if not (vids or attributes):
            vids=True
        if vids:
            header="vid: ["
        else:
            header="[ "
        if self.Type(0)==VariableType.STATE:
            # print the State nature of the 1st variable
            header+=" Optimal State ]"
            comma=False
            if len(smoothed) > 0:
                # Same principle for smoothed probabilities,
                header+="[ "
                for v in smoothed:
                    header+=self.__attributes[v]+" "
                header+="][ "
                filt=lambda var: (var not in smoothed) and (var > 0)
                indices=filter(filt, range(len(self.__attributes)))
                for v in indices:
                    header+=self.__attributes[v]+" "
                header+="]"
            else:
                header+="[ "
        else:
            # separating comma
            header+=self.__attributes[0]
            comma=True
        if len(smoothed) == 0:
            for index in range(len(self.__attributes)-1):
                if comma:
                    header+=", "
                else:
                    comma=True
                header+=self.__attributes[index+1]
            header+=" ]"
##        header is already completed by case 
##        (self.Type(0)==VariableType.STATE) and (len(smoothed) > 0)
##        above
##        else: 
##            # Same principle for smoothed probabilities,
##            header+="[ "
##            for v in smoothed:
##                header+=self.__attributes[v]+" "
##            header+="][ "
##            filt=lambda var: (var not in smoothed) and (var > 0)
##            indices=filter(filt, range(len(self.__attributes)))
##            for v in indices:
##                header+=self.__attributes[v]+" "
##            header+="]"
        return header+"\n"+\
            self.__display(key, vids, attributes, mtg_vids, smoothed, rounded)
    
    def _mtg_tid(self):
        # return the vid in the MTG corresponding to the root of self
        if self.__mtg_tid is None:
            raise Warning, "Current Tree object has not been obtained from " \
                "a MTG"
        else:
            return self.__mtg_tid

    def _mtg_write(self, tabulation="", post_tabulation=""):
        # return a string corresponding to the MTG code for tree self
        # post_tabulation is a string printed between the vertices
        # and the attributes (cf. Trees.Save()) 
        return self.__mtg_skip(self.Root(), tabulation, post_tabulation)

    def _set_int_value_type(self, variable):
        # turns a state variable into an integer-valued variable
        if (variable < 0 or variable >= self.NbVariables()):
            raise IndexError, "variable index out of range: "+str(variable)
        elif (self.__types[variable]==VariableType.STATE):
            self.__types[variable]=VariableType.INT_VALUE
        else:
            msg="bad type for variable " + str(variable) + ": " + \
                str(self.__types[variable])
            raise TypeError, msg
           
    def __display(self, key, vids, attributes, mtg_vids, smoothed, rounded):
        res, lineskip=self.__display_skip(key, "", False, vids, attributes, 
                                          mtg_vids, smoothed, rounded)
        return res

    def __display_skip(self, key, tabulation, lineskip, vids, attributes, 
                       mtg_vids, smoothed, rounded):
        tab= "|-"
        stream=""
        newlineskip=False
        if vids:
            if mtg_vids:
                stream+=str(self.MTGVertex(key))
            else:
                stream+=str(key)
            if attributes:
                stream+=': '
        if attributes:
            if self.Type(0)==VariableType.STATE:
                # A state variable must be written separately
                complete_value=self.Get(key).Values()
                stream+=str([complete_value[0]])
                if len(smoothed) > 0:
                    # Same principle for smoothed probabilities,
                    # using "%.4" % (var)
                    stream+="["
                    for v in smoothed:
                        stream+=str(eval("%.4f" % complete_value[v]))+" "
                    stream+="]"
                    # display other variables except state variable
                    filt1=lambda var: (var not in smoothed) and (var > 0) \
                        and (var not in rounded)
                    filt2=lambda var: (var not in smoothed) and (var > 0) \
                        and (var in rounded)
                    indicesnr=filter(filt1, range(len(complete_value)))
                    indicesr=filter(filt2, range(len(complete_value)))
                    if (len(indicesr) + len(indicesnr)) > 0:
                        stream+='['
                        for v in range(len(complete_value)):
                            if v in indicesnr:
                                stream+=str(complete_value[v]) + " "
                            elif v in indicesr:
                                stream+=\
                                    str(eval("%.4f" % complete_value[v]))+" "
                        stream+=']'
                else:
                    # No smoothed probabilities
                    stream+=str(complete_value[1:len(complete_value)])
            else:
                # No state
                stream+=str(self.Get(key))
            if ((attributes) and not(self.__ctree.IsRoot(key))):
                stream+=self.EdgeType(self.__ctree.Parent(key), key)
        stream+="\n"
        lineskip=False
        children_it= self.__ctree.Children(key)
        children=[] # list of children of vertex key
        current_child=0
        try:
            while 1:
                children.append(children_it.next())
        except StopIteration: pass
        if current_child <= len(children)-1:
        # at least one child remaining
            s=len(tabulation)
            if s > 0:
                tabulation=tabulation[0:len(tabulation)-1]+" "
            tabulation+=tab
            stream+=tabulation
            if current_child == len(children)-1:
            # last child of node v is being written:
            # no vertical bar displayed
                tabulation=tabulation[0:len(tabulation)-2]+"  "
                newlineskip=True
            else:
                newlineskip=False

            restream, lineskip=self.__display_skip(children[current_child],
                                                   tabulation, lineskip,
                                                   vids, attributes, mtg_vids,
                                                   smoothed, rounded)
            stream+=restream
            current_child+=1

            while current_child <= len(children)-1:
                stream+=tabulation
                if current_child == len(children)-1:
                # last child of node v is being written:
                # no vertical bar displayed
                    tabulation=tabulation[0:len(tabulation)-2]+"  "
                    newlineskip=True
                else:
                    newlineskip=False
                restream, lineskip=self.__display_skip(children[current_child],
                                                       tabulation, lineskip,
                                                       vids, attributes, mtg_vids,
                                                       smoothed, rounded)
                stream+=restream
                current_child+=1
            if ((len(tabulation) >= 3) and (not lineskip)):
            # a line is skipped after the last child
            # and tabulation is shifted one step to the left
                lineskip=True
                tabulation=tabulation[0:len(tabulation)-3]+"   "
                stream+=tabulation
                stream+='\n'
        else:
        # current_child > len(children)
            if len(tabulation) >= 2:
                tabulation=tabulation[0:len(tabulation)-2]
        return stream, newlineskip

    def __mtg_header(self, variable_names):
        # print the MTG header
        mtgs='#one-scaled MTG generated by module "trees"\n' 
        mtgs+="CODE:\tFORM-A\n\nCLASSES:\n"
        mtgs+="SYMBOL\tSCALE\tDECOMPOSITION\tINDEXATION\tDEFINITION\n"
        mtgs+="$\t0\tCONNECTED\tFREE      \tIMPLICIT\n"
        mtgs+="V\t1\tNONE     \tFREE      \tEXPLICIT\n\n"
        mtgs+="DESCRIPTION:\n"
        mtgs+="LEFT\tRIGHT\tRELTYPE\tMAX\n"        
        mtgs+="V\tV\t+\t?\n"
        mtgs+="V\tV\t<\t?\n\n"
        mtgs+="FEATURES:\n"
        mtgs+="NAME\tTYPE\n\n"
        for var in range(self.NbVariables()):
            mtgs+=variable_names[var]+"\t"
            if self.Type(var)==VariableType.REAL_VALUE:
                mtgs+="REAL"
            else:
                mtgs+="INT"
            mtgs+="\n"
        mtgs+="\nMTG:\nENTITY-CODE"
        mtgs+="\t"*(self.Depth()+1)
        # this would have to be improved in the case of <
        for var in range(self.NbVariables()):
            mtgs+=variable_names[var]+"\t"
        mtgs+="\n/"
        return mtgs
    
    def __mtg_skip(self, key, tabulation, post_tabulation):
        tab= "\t"
        # position in stream is supposed to be set before function call
        stream="V"+str(key)
        # print the attributes in column self.Depth()+1
        stream+="\t"*(self.__ctree.Order()-self.__ctree.Order(key)+1)
        stream+=post_tabulation
        value=self.Get(key) # attributes
        for var in range(len(value)-1):
            stream+=str(value[var])+"\t"
        stream+=str(value[len(value)-1])
        # print the children
        children_it= self.__ctree.Children(key)
        children=[] # list of children of vertex key
        current_child=0
        try:
            while 1:
                children.append(children_it.next())
        except StopIteration: pass
        # commented lines are useful if we do not want to save the attribute
##        if len(children)>1:
        # tabulation+=tab
##        else:
##            tabulation=""
##        if len(children)>1:        
        stream+="\n"
        while current_child <= len(children)-1:
            etype=self.EdgeType(key, children[current_child])
            if etype=='<':
                tmark="^"
            else:
                tmark=""
                tabulation+=tab
            stream+=tabulation + tmark + etype \
                  +self.__mtg_skip(children[current_child], tabulation, 
                                   post_tabulation)
            current_child+=1
            if etype=='+':
                tabulation=tabulation[0:len(tabulation)-1]
        tabulation=tabulation[0:len(tabulation)-1]


##        while current_child <= len(children)-1:
##            stream+=tabulation+self.EdgeType(key, children[current_child]) \
##                  +self.__mtg_skip(children[current_child], tabulation, 
##                                   post_tabulation)
##            current_child+=1
##        tabulation=tabulation[0:len(tabulation)-1]


##        if len(children)==0:
##            stream+="\n"
        return stream

    def __valid_edge(self, parent, child):
        self.__valid_vid(parent)
        self.__valid_vid(child)
        valid=self.__ctree.IsEdge(parent, child)
        if not(valid):
            msg="pair ("+str(parent)+", "+str(child)+") is not a valid edge"
            raise IndexError, msg

    def __valid_value(self, value):
        tv=TreeValue(value)
        reslist=[]
        msg=""
        status=True
        if len(tv) != self.NbVariables():
            msg+="number of attributes ("+str(self.NbVariables())+") "\
            "is incompatible with argument length ("+str(len(tv))+")"
            raise TypeError, msg
        for var in range(len(tv)):
            if type(tv[var])==int:
                if self.Type(var)==VariableType.REAL_VALUE:
                    reslist.append(tv[var]*1.)
                else:
                    reslist.append(tv[var])
            else:
                if self.Type(var)!=VariableType.REAL_VALUE:
                    msg+="incorrect type for element "+str(var)+" of current "\
                    "argument - type 'int' expected"
                    raise TypeError, msg
                reslist.append(tv[var])
        tv=TreeValue(reslist)
        return TreeValue(tv)
    
    def __valid_vid(self, vid):
        if type(vid) != int:
            msg=str(vid)+" has not the type of a valid vertex identifier."
            raise TypeError, msg
        elif (vid<0 or vid >= self.Size()):
            msg=str(vid)+" is not a valid vertex identifier."
            raise IndexError, msg
        return True
   
    def __str__(self):
        """Display the tree with its attributes.

        Note: the vertex identifiers are not displayed."""
        return self._display(self.Root(), False, True, False)

class TreeStructure:
    """An implementation of tree structures, i.e. trees with no attributes."""

    def __init__(self, arg=None, max_size=I_DEFAULT_TREE_SIZE,
                 max_depth=I_DEFAULT_TREE_DEPTH):
        """Initialize a TreeStructure from a tree object or from a distribution.

        Initialize a TreeStructure object from either:
        - a Tree or a TreeStructure object;
        - a distribution and values for the maximal size and depth."""
        if arg is None:
            self.__tree=ctree.TreeStructure()
        elif issubclass(arg.__class__, Tree):
            # arg is supposed to be a tree...
            ctree_val=arg._ctree()
            self.__tree=ctree.TreeStructure(ctree_val)
        elif issubclass(arg.__class__,TreeStructure):
            #... or a tree structure
            self.__tree=ctree.TreeStructure(arg.__tree)
        else:
            #... or a Distribution
            self.__tree=ctree.TreeStructure(arg, max_size, max_depth)

    def Depth(self):
            """Return the tree depth."""
            return self.__tree.Depth()

    def Display(self, vids=True, edgetypes=True):
        """Display the tree with or without the vertex identifiers (vids)
        and with or without the edge types.

        The vids are displayed if and only if boolean argument vids is True.
        The edge types are displayed if and only if boolean argument edgetypes 
        is True."""
        print self.__display(self.Root(), vids, edgetypes)

    def EdgeType(self, parent, child):
        """Return the type of one given edge (parent, child)."""
        self.__valid_edge(parent, child)
        btype=self.__tree.EdgeType(parent, child)
        if btype:
            return "<"
        else:
            return "+"

    def Root(self):
        """Return the root vertex identifier (vid)."""
        return self.__tree.Root()

    def SelectSubTree(self, vid, keep=True):
        """Select and return a subtree.

        Argument vid, which must be a valid vertex identifier, is the root
        of the selected subtree. Boolean argument keep is True if the
        subtree must be kept, False if it must be pruned from tree."""
        # dummytree=Tree([0], self)
        # res=TreeStructure(dummytree.SelectSubTree(vid, keep))
        # This is an alternative solution. Quite awkward, but
        # does not require select_subtree(Unlabelled_tree,...) to work.

        self.__valid_vid(vid)
        cstr=ctree.SelectSubTreeStructure(self.__tree, vid, keep)
        res=TreeStructure()
        res.__tree=ctree.TreeStructure(cstr)
        return res

    def Simulate(self, distribution, max_size=I_DEFAULT_TREE_SIZE,
                 max_depth=I_DEFAULT_TREE_DEPTH):
        """Simulate a tree structure from a distribution.

        Argument distribution is used for the independant simulation of the
        number of children of each vertex.
        Arguments max_size and max_depth represent the maximal size
        and depth of the simulated tree."""
        self.__tree.Simulate(distribution, max_size, max_depth)

    def Size(self):
        """Return the number of vertices."""
        return self.__tree.Size()

    def _tree(self):
        return self.__tree

    def __display(self, key, vids, edgetypes):
        return self.__display_skip(key, "", False, vids, edgetypes)

    def __display_skip(self, key, tabulation, lineskip, vids, edgetypes):
        tab= "|-"
        stream=""
        if vids:
            stream+=str(key)
        else:
            stream+="*"
        stream+="\n"
        if ((edgetypes) and not(self.__tree.IsRoot(key))):
            stream+=self.EdgeType(self.__tree.Parent(key), key)

        lineskip=False
        children_it= self.__tree.Children(key)
        children=[] # list of children of vertex key
        current_child=0
        try:
            while 1:
                children.append(children_it.next())
        except StopIteration: pass
        if current_child <= len(children)-1:
        # at least one child remaining
            s=len(tabulation)
            if s > 0:
                tabulation=tabulation[0:len(tabulation)-1]+" "
            tabulation+=tab
            stream+=tabulation
            if current_child == len(children)-1:
            # last child of node v is being written:
            # no vertical bar displayed
                tabulation=tabulation[0:len(tabulation)-2]+"  "
            stream+=self.__display_skip(children[current_child],
                                        tabulation, lineskip, vids, edgetypes)
            current_child+=1

            while current_child <= len(children)-1:
                stream+=tabulation
                if current_child == len(children)-1:
                # last child of node v is being written:
                # no vertical bar displayed
                    tabulation=tabulation[0:len(tabulation)-2]+"  "
                stream+=self.__display_skip(children[current_child],
                                            tabulation, lineskip, vids, 
                                            edgetypes)
                current_child+=1
            if ((len(tabulation) >= 3) and (not lineskip)):
            # a line is skipped after the last child
            # and tabulation is shifted one step to the left
                lineskip=True
                tabulation=tabulation[0:len(tabulation)-3]+"   "
                stream+=tabulation
                stream+='\n'
        else:
        # current_child > len(children)
            if len(tabulation) >= 2:
                tabulation=tabulation[0:len(tabulation)-2]
        return stream

    def __valid_vid(self, vid):
        if type(vid) != int:
            msg=str(vid)+" has a not the type of a valid vertex."
            raise TypeError, msg
        elif vid >= self.Size():
            msg=str(vid)+" is not a valid vertex identifier."
            raise IndexError, msg
        return True
    
    def __valid_edge(self, parent, child):
        self.__valid_vid(parent)
        self.__valid_vid(child)
        valid=self.__tree.IsEdge(parent, child)
        if not(valid):
            msg="pair ("+str(parent)+", "+str(child)+") is not a valid edge"
            raise IndexError, msg
        else:
            return True
        
    def __valid_value(self, value):
        tv=TreeValue(value)
        reslist=[]
        msg=""
        status=True
        if len(tv) != self.NbVariables():
            msg+="number of attributes ("+str(self.NbVariables())+") "\
            "is incompatible with argument length ("+str(len(tv))+")"
            raise TypeError, msg
        for var in range(len(tv)):
            if type(tv[var])==int:
                if self.Type(var)==VariableType.REAL_VALUE:
                    reslist.append(tv[var]*1.)
                else:
                    reslist.append(tv[var])
            else:
                if self.Type(var)!=VariableType.REAL_VALUE:
                    msg+="incorrect type for element "+str(var)+" of current "\
                    "argument - type 'int' expected"
                    raise TypeError, msg
                reslist.append(tv[var])
        tv=TreeValue(reslist)
        return TreeValue(tv)

    def __str__(self):
        return str(self.__tree)



class Trees(object):
    """A set of trees with integral and floating attributes."""

          
    def __init__(self, arg, arg2=None, 
                 attribute_names=None, attribute_def=None, scale=None):
        """Initialize a Trees object from trees.

        Initialize a Trees object from either:
        - a list of Tree objects;
        - a Trees object.
        - a MTG file, a filter on the vertices, a list of attribute names,
          a list of attribute functions and the considered scale."""
        self.__mtg_to_tree_vid=None
        self.__tree_to_mtg_vid=None
        self.__tree_to_mtg_tid=None
        self.__mtg_to_tree_tid=None
        if issubclass(arg.__class__, Trees):
            # arg is supposed to be a Trees object...
            self.__ctrees=ctrees.CTrees(arg.__ctrees)
            self.__types=list(arg.__types)
            self.__tmap=list(arg.__tmap)
            self.__tmapi=list(arg.__tmapi)
            self.__attributes=list(arg.__attributes)
            arg._copy_vid_conversion(self)
            arg._copy_tid_conversion(self)
        elif issubclass(arg.__class__, ctrees.CTrees):
            # ... or a CTrees object...
            self.__ctrees=ctrees.CTrees(arg)
            if arg2 is None:
                # arg2 corresponds to the types in this case
                nbvariables=arg.NbInt()+arg.NbFloat()
                self.__types=[]
                for var in range(nbvariables):
                    self.__types.append(VariableTypeDict[arg.Type(var)])
            else:
                for var in range(len(arg2)):
                # checking that the type is valid
                    if not (arg2[var]==VariableType.INT_VALUE
                        or arg2[var]==VariableType.STATE
                        or arg2[var]==VariableType.REAL_VALUE):
                        msg="unknown type : "+str(arg2[var])
                        raise TypeError, msg
                self.__types=list(arg2)
        elif type(arg)==str:
            # ... or the name of a MTG file...
            self.__mtg_to_tree_vid=[]
            # conversion from MTG vid to tree vid
            self.__tree_to_mtg_vid=[]
            # conversion from tree vid to MTG vid
            self.__tree_to_mtg_tid={}
            # conversion from tree id to MTG ComponentRoots
            self.__mtg_to_tree_tid={}
            # conversion from MTG ComponentRoots to tree id
            nbmtg=0     # number of MTG ComponentRoots
            nbtrees=0   # number of trees
            no_scale=(scale is None)
            import openalea.aml as amlPy
            mode=False
            if not amlPy.getmode():
            # conversion from AML object to Python
                mode=True
                amlPy.setmode(1)
            try:
                mtg_file=file(arg, 'r')
            except IOError, msg:
                raise IOError, msg
            M=amlPy.MTG(arg)
            # reading the MTG header
            nbfloat=0
            nbint=0
            mtg_code=mtg_file.readline()
            # seeking for keyword "SCALE"
            found=False
            while not found:
                mtg_code=mtg_file.readline()
                pos=mtg_code.upper().find("SCALE")
                if pos != -1:
                    # keyword "SCALE" has been found but can be commented out
                    if (mtg_code.find("#",0,pos)==-1):
                        found=True
            mtg_code=mtg_file.readline()
            while mtg_code.isspace():
                mtg_code=mtg_file.readline()
            while ((mtg_code.upper().find("DESCRIPTION")==-1) and
                    not mtg_code.isspace()):
                words=mtg_code.split()
                if no_scale:
                # scale should not be modified if given by user
                    try:
                        scale=int(words[1])
                    except ValueError: pass
                mtg_code=mtg_file.readline()
            # same procedure for keyword "FEATURES"
            found=False
            while not found:
                mtg_code=mtg_file.readline()
                pos=mtg_code.upper().find("FEATURES")
                if pos != -1:
                    # keyword "FEATURES" has been found but can be commented out
                    if (mtg_code.find("#",0,pos)==-1):
                        found=True
            # same procedure for keyword "NAME"
            found=False
            while not found:
                mtg_code=mtg_file.readline()
                pos=mtg_code.upper().find("NAME")
                if pos != -1:
                    # keyword "NAME" has been found but can be commented out
                    if (mtg_code.find("#",0,pos)==-1):
                        # keyword "NAME" must be followed by TYPE
                        if (mtg_code.find("TYPE",pos)!=-1):
                            found=True
            mtg_code=mtg_file.readline()
            while mtg_code.isspace():
                mtg_code=mtg_file.readline()
            feature_names=[]
            feature_types=[]
            self.__types=[]
            # determining the type of each feature 
            # if attribute_def is unspecified
            if attribute_def is None:
                while ((mtg_code.upper().find("MTG")==-1) and
                        not mtg_code.isspace()):
                    words=mtg_code.split()
                    if words[1].upper()=='INT':
                        feature_names.append(words[0])
                        feature_types.append(VariableType.INT_VALUE)
                        nbint+=1
                    elif words[1].upper()=='REAL':
                        feature_names.append(words[0])
                        feature_types.append(VariableType.REAL_VALUE)
                        nbfloat+=1
                    mtg_code=mtg_file.readline()
            mtg_file.close()
            if attribute_def is None:
                if attribute_names is None:
                    attribute_names=feature_names
                attribute_def=[]
                for name in attribute_names:
                    attribute_def.append(lambda x,y=name: amlPy.Feature(x,y))
                self.__types=feature_types
            elif attribute_names is None:
                # only attribute_def is given
                attribute_names=[]
                for var in range(len(attribute_def)):
                    attribute_names.append("Variable"+str(var))
            else:
                # both attribute_def and attribute_name are given
                index=0
                if len(attribute_names)!=len(attribute_def):
                    msg="inconsistent number of attributes in name list and" \
                        +" definition list"
                    raise ValueError, msg
                for name in attribute_names:
                    index+=1
                    if (type(name)!=str):
                        msg="bad type for name of attribute "+str(index) \
                            + ": type 'str' expected"
                        raise TypeError, msg
                index=0
                for func in attribute_def:
                    if not callable(func):
                        msg="bad type for function defining attribute "\
                            +str(index)+": expecting a callable"
                        raise TypeError, msg
                    index+=1                    
            if arg2 is None:
            # default filter: accept all vertices
                arg2=lambda x: True
            elif not callable(arg2):
                msg="bad type for filtering function: expecting a callable"
                raise TypeError, msg
            croots=amlPy.ComponentRoots(amlPy.MTGRoot(), Scale=scale)
            v0=croots[0]
            tree_list=[]
            for v in croots:
                # roots at given scale. Each value of v represents a new tree.
                nbmtg+=1
                if ((len(self.__types)==0) and (v0==v)):
                    # determine the type of each variable if necessary
                    values=[attr(v) for attr in attribute_def]
                    index=0
                    for val in values:
                        if type(val)==int:
                            self.__types.append(VariableType.INT_VALUE)
                            nbint+=1
                        elif type(val)==float:
                            self.__types.append(VariableType.REAL_VALUE)
                            nbfloat+=1
                        else:
                            msg="bad result type for function defining attribute "\
                                +str(index)+": expecting types 'int' or 'float'"
                            raise TypeError, msg
                        index+=1
                if (v0==v):
                    # build a list of typical values
                    typical_values=[]
                    for t in self.__types:
                        if ((t==VariableType.INT_VALUE) or 
                            (t==VariableType.STATE)):
                            typical_values.append(0)
                        else:
                            typical_values.append(0.)
                mtg2tree_vertices={}
                # dictionnary vid(mtg)->vid(Tree)
                mtgvertices=amlPy.Descendants(v) # warning : contains v !
                current_tree=ctree.CTree(nbint, nbfloat, 0, 0)
                for x in mtgvertices:
                    if arg2(x):
                        values=[attr(x) for attr in attribute_def]
                        typecheck=[((type(val)==int) or (type(val)==float))
                                   for val in values]
                        if typecheck.count(0) > 0:
                            # at least one element of values has incorrect type
                            msg="bad type for attributes of vertex "\
                                +str(x)+": expecting types 'int' or 'float'"
                            raise TypeError, msg                     
                        mtg2tree_vertices[x]= \
                            current_tree.AddVertex(ctree.TreeValue(values))
                        if ((x != v) and mtg2tree_vertices.has_key(amlPy.Father(x))):
                            edgetype=amlPy.EdgeType(amlPy.Father(x), x)
                            if edgetype=='+':
                                btype=False
                            elif edgetype=='<':
                                btype=True
                            elif type(edgetype)==str:
                                msg="Bad value for edge type: " + edgetype
                                raise ValueError, msg
                            else:
                                msg="Bad type for edge type: " + str(type(edgetype))
                                raise TypeError, msg
                            current_tree.AddEdge(mtg2tree_vertices[amlPy.Father(x)],
                                                 mtg2tree_vertices[x], btype)
                        elif (x != v):
                            # the father of x has been filtered but x has not
                            msg="father of vertex "+str(x)+" has been deleted by " \
                                "filter but vertex "+str(x)+" has not"
                            raise IndexError, msg
                if current_tree.Size() > 0:
                    tree_list.append(Tree(typical_values, current_tree))
                    nbtrees+=1
                    self.__tree_to_mtg_tid[nbtrees-1]=v # nbmtg-1
                    # conversion of the MTG vid to Trees vid and vice versa
                    self.__mtg_to_tree_vid.append(mtg2tree_vertices)
                    tree2mtg_vertices={}
                    for key in mtg2tree_vertices.keys():
                        tree2mtg_vertices[mtg2tree_vertices[key]]=key
                    self.__tree_to_mtg_vid.append(tree2mtg_vertices)            
            if len(tree_list) > 0:
                self.__ctrees=ctrees.CTrees(tree_list)
            else:
                raise ValueError, "cannot build an empty Trees object"
            # inversion of table of the tree_to_mtg identifiers
            for key in self.__tree_to_mtg_tid.keys():
                self.__mtg_to_tree_tid[self.__tree_to_mtg_tid[key]]=key
            if len(self.__types)==0:
                for val in values:
                    if type(val)==float:
                        self.__types.append(VariableType.INT_VALUE)
                    else:
                        self.__types.append(VariableType.REAL_VALUE)
            if mode:
                amlPy.setmode(0)
        else:
            #... or a list of Tree objects
            self.__ctrees=ctrees.CTrees(arg)
            self.__types=list(arg[0].Types())
            if attribute_names is None:
                attribute_names=list(arg[0].Attributes())
            build_mtg_to_tree=False
            for t in range(len(arg)):
                try:
                    arg[t].TreeVertex()
                except Warning:
                    pass
                else:
                    build_mtg_to_tree=True
                    break
            if build_mtg_to_tree:
                self.__mtg_to_tree_vid=[]
                self.__mtg_to_tree_tid={}
                for t in range(len(arg)):
                    try:
                        tv=arg[t].TreeVertex()
                        tid=arg[t]._mtg_tid()
                    except Warning:
                        self.__mtg_to_tree_vid.append({})
                    else:
                        self.__mtg_to_tree_vid.append(dict(tv))
                        self.__mtg_to_tree_tid[tid]=t

            build_tree_to_mtg=False
            for t in range(len(arg)):
                try:
                    arg[t].MTGVertex()
                except Warning:
                    pass
                else:
                    build_tree_to_mtg=True
                    break
            if build_tree_to_mtg:
                self.__tree_to_mtg_vid=[]
                self.__tree_to_mtg_tid={}
                for t in range(len(arg)):
                    try:
                        tv=arg[t].MTGVertex()
                        tid=arg[t]._mtg_tid()
                    except Warning:
                        self.__tree_to_mtg_vid.append({})
                    else:
                        self.__tree_to_mtg_vid.append(dict(tv))
                        self.__tree_to_mtg_tid[t]=tid
                
            # build the dictionnaries corresponding to the tree -> MTG
            # and MTG -> tree vid conversions

        if not issubclass(arg.__class__, Trees):
            # The variables in CTrees are packed as follows:
            # [integral_variables, real_variables].
            # Thus, self.__tmap maps the variables of self to
            # the variables of self.__ctrees
            self.__tmap=range(len(self.__types))
            nbint=0
            for var in range(len(self.__types)):
                if self.__types[var]!=VariableType.REAL_VALUE:
                    nbint+=1
            findex=nbint
            iindex=0
            for var in range(len(self.__types)):
                if self.__types[var]!=VariableType.REAL_VALUE:
                    self.__tmap[var]=iindex
                    iindex+=1
                else:
                    self.__tmap[var]=findex
                    findex+=1
            # computation of the inverse permutation of self.__tmap
            self.__tmapi=range(len(self.__tmap))
            for var in range(len(self.__tmap)):
                self.__tmapi[self.__tmap[var]]=var
        if attribute_names is None:
            attribute_names=[]
            for var in range(self.NbVariables()):
                attribute_names.append("Variable"+str(var))
        self.__attributes=list(attribute_names)
        if len(self.__attributes) != len(self.__types):
            msg="Number of variables ("+str(len(self.__types))+\
                ") and number of attribute names ("+str(len(self.__attributes))+\
                ") do not match"
            raise Warning, msg

    def Attributes(self):
        """Return the name of the tree attributes."""
        return list(self.__attributes)

    def BuildSequences(self, maximal_sequences=True):
        """Extract Sequences from the Trees, 
        cutting or not sequences after branching"""
        import os
        # print the sequences into a file
        prefix="seqtmp"
        file_created=False
        while not file_created:
            try:
                cfile=open(prefix+'.seq','r')
            except IOError:
                # file does not exist
                file_name= prefix+".seq"
                file_created=True
            else:
                cfile.close()
                import random
                prefix+=str(random.randint(1,9))                
        try:
            self.__ctrees.BuildSequences(file_name, maximal_sequences)
        except RuntimeError, error:
            os.remove(file_name)
            raise FormatError(error)
        else:
            import openalea.aml as amlPy
            res= amlPy.Sequences(file_name)
            os.remove(file_name)            
            return res

    def BuildPySequences(self, maximal_sequences=True):
        """Extract sequence_analysis.Sequences from the Trees,
        cutting or not sequences after branching"""
        import os
        try:
            res = self.__ctrees.BuildPySequences(maximal_sequences)
        except RuntimeError, error:
            raise FormatError(error)
        else:
            return res.markovian_sequences()

    def BuildVectors(self):
        """Extract Vectors from the Trees."""
        import os
        return self.__ctrees.BuildVectors()

    def Cluster(self, mode, variable, limit):
        """Clustering of values.

        Cluster the values of a given variable.
        Argument mode must be "Step" or "Limit" and corresponds to the
        clustering mode.
        Argument variable refers to the clustered variable.
        Argument limit can be either:
        - the clustering step if mode is "Step"
        - the list of bounds defining the clusters if mode is "Limit". """
        cvariable=self._valid_cvariable(variable)+1
        # correspondence of variables between self and self.__ctrees
        try:
            if not (mode.upper()=="STEP" or mode.upper()=="LIMIT"):
                msg="unknown clustering mode: "+str(mode)\
                    +" - expecting 'Step' or 'Limit'"
                raise ValueError, msg
            elif (mode.upper()=="STEP" and type(limit)!=int):
                msg="incorrect type for clustering step: "+str(type(limit))
                raise TypeError, msg
            elif (mode.upper()=="LIMIT" and 
                  not (hasattr(limit, "__getitem__"))):
                msg="incorrect type for clustering limits: "+str(type(limit))
                raise TypeError, msg                                
            else:
                cclustered=self.__ctrees.Cluster(cvariable, limit)
        except RuntimeError, error:
            # raise stat_tool.Format_error, str(cerror)
            replaced=str(error) # self.__replacestr(str(cerror), "variable")
            raise FormatError, replaced
        clustered=Trees(cclustered, self.__types, self.__attributes)
        self._copy_vid_conversion(clustered)
        self._copy_tid_conversion(clustered)
        return clustered

    def ComputeStateTrees(self, model, algorithm="Viterbi", characteristics=True):
        """Compute the optimal state trees corresponding to the observed trees.

        :Parameters:

          * `model` (`hiddenMarkovtree`) - model used for the computation,
          * `algorithm` (str) - type of algorithm ("Viterbi" or "ForwardBackward"),
          * `characteristics` (bool) - characteristic distributions are computed iif \
            argument is True.
        """
        import openalea.tree_statistic.hmt, openalea.tree_statistic.hmt.chmt
        hmt=openalea.tree_statistic.hmt
        chmt=openalea.tree_statistic.hmt.chmt
        RestorationAlgorithm=stat_tool.RestorationAlgorithm
        if not issubclass(model.__class__, hmt.HiddenMarkovTree):
            msg='bad type for argument "hmt": HiddenMarkovTree expected'
            raise TypeError, msg
        if type(algorithm)!=str:
            msg='bad type for argument "algorithm:"'+"type 'str' expected"
            raise TypeError, msg
        elif ((algorithm.upper()!="FORWARDBACKWARD")
              and (algorithm.upper()!="VITERBI")):
            msg='bad value for argument "algorithm":'+arg2
            raise ValueError, msg
        if (algorithm.upper()=="FORWARDBACKWARD"):
            arg5="FORWARD_BACKWARD"                    
        algorithm="RestorationAlgorithm."+algorithm.upper()
        StateTrees=eval(algorithm)
        chmt_data=(model._chmt()).ComputeStateTrees(self.__ctrees, StateTrees,
                                                    characteristics)
        chmt_data=chmt_data.StateTrees()
##        smoothed_names=["Smoothed"+str(s) 
##                        for s in range(model._chmt().NbStates())]
        res=hmt.HiddenMarkovTreeData(chmt_data, model._chmt(), aliasing=True, 
                                     attribute_names=["OptimalState"]
                                     +self.__attributes) #+smoothed_names)
        self._copy_vid_conversion(res)
        self._copy_tid_conversion(res)
        return res

    def Difference(self, variable=None):
        """First-order differentiation of trees."""
        if (variable is None):
            cvariable=stat_tool.I_DEFAULT()
        else:
            cvariable=self._valid_cvariable(variable)+1
        try:
            cdiff=self.__ctrees.Difference(cvariable)
        except RuntimeError, error:
            replaced=self.__replacestr(str(error), "variable")
            raise FormatError, replaced
        diff=Trees(cdiff, self.__types, self.__attributes)
        self._copy_vid_conversion(diff)
        self._copy_tid_conversion(diff)
        return diff
            
    def Display(self, ViewPoint=None, Detail=1):
        """Display Trees object with a level of Detail 1 or 2.
        
        Usage: Display(Detail=2)
               Display(ViewPoint="Data", Detail=1)"""
        if Detail==1:
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
            
        replaced=self.__replacestr(self.__ctrees.Display(exhaustive), 
                                   "variable", True)
        print replaced
        if not (ViewPoint is None):
            if type(ViewPoint)!=str:
                msg="Bad type for 'ViewPoint' argument:"+str(type(Detail)) \
                +" - expecting type 'str'"
                raise TypeError, msg
            elif ViewPoint.upper()=="DATA":
                for t in range(self.NbTrees()):
                    print "("+str(t)+")"
                    self.Tree(t).Display(vids=True, attributes=True)
            else:
                msg="Bad value for 'ViewPoint' argument:"+str(type(Detail)) \
                +" - expecting 'DATA'"
                raise ValueError, msg
        
    def Estimate(self, model_name, arg1, arg2=None, arg3=None, arg4=None, 
                 arg5=None, arg6=None, Algorithm="ForwardBackward", Saem=1., 
                 ForceParametric=[]):
        """Estimate a (hidden) Markov tree.
        
        Algorithm correspond to the type of restoration/maximisation algorithm:
        'ForwardBackward', 'Viterbi', 'ForwardBackwardSampling'
        or 'GibbsSampling'
        Saem correspond to the rate of decay of the part corresponding 
        to restored states. Saem=0. for pure SEM or CEM algorithms.
        
        :Usage:

            Estimate("HIDDEN_MARKOV_TREE", nb_state, structure,
                          InitialSelfTransition, NbIteration, StateTrees, 
                          Algorithm, Saem, Counting, ForceParametric)
            Estimate("HIDDEN_MARKOV_TREE", hmt, NbIteration, Algorithm, 
                          Saem, Counting)"""
        import openalea.tree_statistic.hmt, openalea.tree_statistic.hmt.chmt
        hmt=openalea.tree_statistic.hmt
        chmt=openalea.tree_statistic.hmt.chmt
        RestorationAlgorithm=stat_tool.RestorationAlgorithm
        chmt_data=chmt.CHmt_data(self._ctrees())
        if type(model_name)==str:
            if type(Algorithm)!=str:
                msg='bad type for argument "Algorithm" in ' \
                'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                'SelfTransition, NbIteration, StateTrees, Algorithm, '\
                'Saem, Counting): '\
                +"type 'str' expected"
                raise TypeError, msg                
            if Algorithm.upper()=="FORWARDBACKWARD":
                algo="FORWARD_BACKWARD"
                algo="RestorationAlgorithm."+algo.upper()
            elif Algorithm.upper()=="FORWARDBACKWARDSAMPLING":
                algo="FORWARD_BACKWARD_SAMPLING"
                algo="RestorationAlgorithm."+algo.upper()
            elif Algorithm.upper()=="GIBBSSAMPLING":
                algo="GIBBS_SAMPLING"
                algo="RestorationAlgorithm."+algo.upper()
            elif Algorithm.upper()=="VITERBI":
                algo="VITERBI"
                algo="RestorationAlgorithm."+algo.upper()
            else:
                msg='bad value for argument "Algorithm" in ' \
                'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                  'SelfTransition, NbIteration, StateTrees, Algorithm, '\
                  'Saem, Counting): '\
                  +str(Algorithm)
                raise ValueError, msg
            EMAlgo=eval(algo)
            if type(Saem)!=float:
                msg='bad type for argument "Saem" in ' \
                'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                'SelfTransition, NbIteration, StateTrees, Algorithm, '\
                'Saem, Counting): '\
                +"type 'float' expected"
                raise TypeError, msg
            if model_name.upper()=="HIDDEN_MARKOV_TREE":
                if type(arg1)==int:
                    # Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, 
                    #          InitialSelfTransition, NbIteration, StateTrees, 
                    #          Counting, Algorithm, Saem, ForceParametric)
                    if arg2 is None:
                        msg='argument "structure" is mandatory in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric)'
                        raise TypeError, msg
                    elif type(arg2)!=str:
                        msg='bad type for argument "structure" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +"type 'str' expected"
                        raise TypeError, msg
                    elif ((arg2.upper()!="IRREDUCTIBLE")
                          and (arg2.upper()!="IRREDUCIBLE")
                          and (arg2.upper()!="LEFTRIGHT")):
                        msg='bad value for argument "structure" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +arg2+" - expecting 'Irreducible' or 'LeftRight'"
                        raise ValueError, msg
                    structure=(arg2.upper()=="LEFTRIGHT")
                    if arg3 is None:
                        arg3=stat_tool.SELF_TRANSITION
                    elif type(arg3)!=float:
                        msg='bad type for argument "SelfTransition" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +"type 'float' expected"
                        raise TypeError, msg                        
                    if arg4 is None:
                        arg4=stat_tool.D_DEFAULT
##                        msg='argument "NbIteration" is mandatory in ' \
##                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
##                        'SelfTransition, NbIteration, StateTrees, Counting, '\
##                        'Algorithm, Saem, ForceParametric): '
##                        raise TypeError, msg
                    elif type(arg4)!=int:
                        msg='bad type for argument "NbIteration" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +"type 'int' expected"
                        raise TypeError, msg
                    if arg5 is None:
                        arg5="VITERBI"
                    elif type(arg5)!=str:
                        msg='bad type for argument "StateTrees" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +"type 'str' expected"
                        raise TypeError, msg
                    elif ((arg5.upper()!="FORWARDBACKWARD")
                          and (arg5.upper()!="VITERBI")):
                        msg = 'bad type for argument "StateTrees" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '\
                        +arg2
                        raise ValueError, msg
                    if (arg5.upper()=="FORWARDBACKWARD"):
                        arg5="FORWARD_BACKWARD"
                    arg5="RestorationAlgorithm."+arg5.upper()
                    StateTrees=eval(arg5)
                    if arg6 is None:
                        arg6=True
                    elif type(arg6)!=int:
                        msg='bad type for argument "Counting" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", nb_state, structure, '\
                        'SelfTransition, NbIteration, StateTrees, Counting, '\
                        'Algorithm, Saem, ForceParametric): '
                        +"boolean type expected"
                        raise TypeError, msg
                    try:
                        chmt=chmt_data.EstimationCiHmot(arg1, structure, arg6, 
                                                        StateTrees, EMAlgo, 
                                                        Saem, arg3, arg4, 
                                                        ForceParametric)
                    except RuntimeError, error:
                        raise FormatError(error)
                elif issubclass(arg1.__class__, hmt.HiddenMarkovTree):
                    # Estimate("HIDDEN_MARKOV_TREE", hmt, NbIteration, Counting,
                    #          Algorithm, Saem, ForceParametric)
                    if arg2 is None:
                        arg2=stat_tool.D_DEFAULT
##                        msg='argument "NbIteration" is mandatory in ' \
##                        'Estimate("HIDDEN_MARKOV_TREE", hmt, ' \
##                                  'NbIteration, Algorithm, Saem, Counting)'
##                        raise TypeError, msg
                    elif type(arg2)!=int:
                        msg='bad type for argument "NbIteration" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", hmt, NbIteration, ' \
                                  'Counting, Algorithm, Saem, ForceParametric): '\
                          +"type 'int' expected"
                        raise TypeError, msg
                    if arg3 is None:
                        arg3=True
                    elif ((type(arg3)!=int) and (type(arg3)!=bool)):
                        msg='bad type for argument "Counting" in ' \
                        'Estimate("HIDDEN_MARKOV_TREE", hmt, NbIteration, ' \
                                  'Counting, Algorithm, Saem, ForceParametric): '\
                          +"boolean type expected"
                        raise TypeError, msg
                    if (ForceParametric==[]):
                        ForceParametric=False
                    elif type(ForceParametric)!=bool:
                        ForceParametric=False
                    try:
                        chmt=chmt_data.EstimationCiHmot(arg1._chmt(), arg3, 
                                                        RestorationAlgorithm.VITERBI,
                                                        EMAlgo, Saem, arg2, 
                                                        ForceParametric)
                    except RuntimeError, error:
                        raise FormatError(error)
                else:
                    msg="bad type for argument 1: "+type(arg1)
                    raise TypeError, msg
            else:
                msg="unknown model name: "+model_name
                raise ValueError, msg
        else:
            raise TypeError, "bad type for argument 1: type 'str' expected"
        estimated_hmt=hmt.HiddenMarkovTree(chmt, True)
        self._copy_vid_conversion(estimated_hmt)
        self._copy_tid_conversion(estimated_hmt)
        estimated_hmt._attributes=self.Attributes()
        return estimated_hmt
        
    def ExtractHistogram(self, nature, variable=None, value=None):
        """Extract a frequency distribution from the Trees.
        
        Usage:  ExtractHistogram("Size")
                ExtractHistogram("NbChildren")
                ExtractHistogram("Value", variable)
                ExtractHistogram("VariableName")
                ExtractHistogram("NbZones", variable, value)"""
        if type(nature)==str:
            if (nature.upper()=="SIZE"):
                chisto=self.__ctrees.ExtractSizeHistogram()
                return chisto
            elif (nature.upper()=="NBCHILDREN"):
                chisto=self.__ctrees.ExtractNbChildrenHistogram()
                return chisto
            elif (nature.upper()=="VALUE"):
                if variable is None:
                    if self.NbVariables()==1:
                    # variable argument is not mandatory if there is only
                    # one variable
                        variable=0
                    else:
                        raise TypeError, 'argument 2 is mandatory in ' \
                                         'ExtractHistogram("Value", variable)'
                try:
                    chisto=self.__ctrees.ExtractValueHistogram(
                        self._valid_cvariable(variable)+1)
                except RuntimeError, error:
                    raise FormatError(error)
                else:
                    return chisto
            elif ((nature.upper()=="FIRSTOCCURRENCEROOT") or
                  (nature.upper()=="FIRSTOCCURRENCELEAVES") or
                  (nature.upper()=="SOJOURNSIZE") or
                  (nature.upper()=="NBZONES") or
                  (nature.upper()=="NBOCCURRENCES")):
                if variable is None:
                    if self.NbVariables()==1:
                    # variable argument is not mandatory if there is only
                    # one variable
                        variable=0
                    else:
                        msg='argument 2 is mandatory in ExtractHistogram' + \
                            '(FeatureName, variable, value)'
                        raise TypeError, msg
                if value is None:
                    # value argument is mandatory 
                    raise TypeError, 'argument 3 is mandatory in ' \
                        'ExtractHistogram(FeatureName, variable, value)'
                elif type(value)!=int:
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
                chartype+=0
                try:
                    chisto=self.__ctrees.ExtractFeatureHistogram(
                        chartype, self._valid_cvariable(variable)+1, value)
                except RuntimeError, error:
                    raise FormatError(error)
                else:
                    return chisto
            else:
                try:
                    v=self.__name_to_variable(nature)
                except ValueError:
                    msg="unknown feature histogram: "+nature
                    raise ValueError, msg
                chisto=self.__ctrees.ExtractValueHistogram(v+1)
                return chisto
        else:
            raise TypeError, "bad type for argument 1: type 'str' expected"

    def Merge(self, tree_list):
        """Merge Trees objects (contained in a list) with self.

        If the argument is a list of Trees objects and if the variables of
        each Trees object are compatible, return a Trees object."""
        #cerror=stat_tool.Format_error()
        ctree_list=[]
        for t in range(len(tree_list)):
            if issubclass(tree_list[t].__class__, Trees):
                ctree_list.append(tree_list[t].__ctrees)
            else:
                ctree_list.append(tree_list[t])
        try:
            # cmerged=self.__ctrees.Merge(cerror, ctree_list)
            cmerged=self.__ctrees.Merge(ctree_list)
        except RuntimeError, error:
            # raise stat_tool.Format_error, str(cerror)
            replaced=self.__replacestr(str(error), "variable")
            # replaced=self.__replacestr(str(cerror), "variable")
            raise FormatError, replaced
            # raise FormatError(replaced)
        merged=Trees(cmerged, self.__types, self.__attributes)
        # self._copy_vid_conversion(merged)
        # self._copy_tid_conversion(merged)
        return merged

    def MPlot(self, ViewPoint="Data", Length=None, BottomDiameter=None, 
              Color=None, DressingFile=None, Title="", variable=0):
        """Graphical output using the Geom 3D viewer for Trees 
           or MultiPlotSet for features.
        
        Usage:  MPlot(ViewPoint="Data")
                MPlot("FirstOccurrenceRoot", variable=0)
        Other possible values for ViewPoint: 
            FirstOccurrenceLeaves
            SojournSize
            Counting"""
        if type(ViewPoint)!=str:
            msg='bad type for argument "ViewPoint": '\
              +"type 'str' expected"
            raise TypeError, msg
        import os
        if ViewPoint.upper()=="DATA":
            # Graphical output of the Trees using the Geom 3D viewer 
            # create the MTG file
            mtgprefix="gvtmp"
            file_created=False
            while not file_created:
                try:
                    cfile=open(mtgprefix+'.mtg','r')
                except IOError:
                    # file does not exist
                    mtgfile_name= mtgprefix+".mtg"
                    file_created=True
                else:
                    cfile.close()
                    import random
                    mtgprefix+=str(random.randint(1,9))
    
            if DressingFile is None:
                # create the dressing file
                drfprefix="dftmp"
                file_created=False
                while not file_created:
                    try:
                        cfile=open(drfprefix+'.drf','r')
                    except IOError:
                        # file does not exist
                        drffile_name= drfprefix+".drf"
                        file_created=True
                    else:
                        cfile.close()
                        import random
                        drfprefix+=str(random.randint(1,9))
                
                dressing=file(drffile_name,'w+')
                dressing.write("Phyllotaxy = 103\n")
                dressing.write("NbPlantsPerLine = "+str(self.NbTrees())+"\n")
                dressing.close()
            else:
                # use the given dressing file
                drffile_name=DressingFile
    
            # create the temporary MTG file
            self.Save(mtgfile_name, False, list(self.__attributes))
    
            import openalea.aml as amlPy
            
            mode=False
            if not amlPy.getmode():
            # conversion from AML object to Python
                mode=True
                amlPy.setmode(1)
                
            M=amlPy.MTG(mtgfile_name)
            DR=amlPy.DressingData(drffile_name)
                            
            if Length is None:
                # define a default length function
                default_lengthfunc= lambda x: 20
                for var in range(self.NbVariables()):
                    lvariable_name=self.__attributes[var]
                    if ((lvariable_name.upper()).find("LEN") != -1):
                        # current variable name contains "LEN"
                        default_lengthfunc= \
                            lambda x: amlPy.Feature(x, lvariable_name)
                        break
            else:
                default_lengthfunc= Length
            # check the type and value
            def lengthfunc(x):
                res=default_lengthfunc(x)
                if (type(res)==int) or (type(res)==float):
                    if (res > 0):
                        return res
                    else:
                        return 10
                else:
                    return 10
    
            if BottomDiameter is None:
                # define a default diameter function
                default_diamfunc= lambda x: 5
                for var in range(self.NbVariables()):
                    dvariable_name=self.__attributes[var]
                    if ((dvariable_name.upper()).find("DIA") != -1):
                        # current variable name contains "DIA"
                        default_diamfunc= \
                            lambda x: amlPy.Feature(x, dvariable_name)
                        break
            else:
                default_diamfunc= BottomDiameter
            # check the type and value
            def diamfunc(x):
                res=default_diamfunc(x)
                if (type(res)==int) or (type(res)==float):
                    if (res > 0):
                        return res
                    else:
                        return 10
                else:
                    return 10
    
            if Color is None:
                # define a default color function
                default_colorfunc= lambda x: 0
                for var in range(self.NbVariables()):
                    cvariable_type=self.__types[var]
                    cvariable_name=self.__attributes[var]
                    if (cvariable_type == VariableType.STATE):
                        # current variable is a state variable
                        default_colorfunc= \
                            lambda x: amlPy.Feature(x, cvariable_name)+2
                        break
            else:
                default_colorfunc= Color
            # check the type and value
            def colorfunc(x):
                res=default_colorfunc(x)
                if (type(res)==int):
                    if (res > 0):
                        return res
                    else:
                        return 0
                else:
                    return 0
    

            vtx_list=amlPy.VtxList(Scale=1)
            pf=amlPy.PlantFrame(vtx_list, Scale=2, DressingData=DR, 
                                Length=lengthfunc, BottomDiameter=diamfunc)
            amlPy.Plot(pf, Color=colorfunc)
            
            if mode:
                amlPy.setmode(0)
            
            # remove the temporary files
            os.remove(mtgfile_name)
            if DressingFile is None:
                os.remove(drffile_name)
                
        else:
            # Graphical output of the features using Gnuplot.py
            if ViewPoint.upper()=="FIRSTOCCURRENCEROOT":
                ftype=CharacteristicType.FIRST_OCCURRENCE_ROOT
            elif ViewPoint.upper()=="FIRSTOCCURRENCELEAVES":
                ftype=CharacteristicType.FIRST_OCCURRENCE_LEAVES
            elif ViewPoint.upper()=="SOJOURNSIZE":
                ftype=CharacteristicType.SOJOURN_SIZE
            elif ViewPoint.upper()=="COUNTING":
                ftype=CharacteristicType.NB_ZONES
            else:
                msg='bad value for argument "ViewPoint": '+ViewPoint
                raise ValueError, msg
            cvariable = self._valid_cvariable(variable)
            print "Using cvariable", cvariable
            if (not self.__ctrees.IsCharacteristic(cvariable, ftype)):
                msg = "Characteristic " + ViewPoint + " not computed " + \
                    "for variable ", variable
                raise FormatError(msg)
            file_id = str(self._valid_cvariable(variable)+1)+str(ftype+1)
            try:
                self.__ctrees.plot(Title=Title, Suffix=file_id,
                                   Params=(ftype, variable))
            except RuntimeError, f:
                raise FormatError(f)

    def MergeVariable(self, tree_list):
        """Merge the variables of Trees objects (contained in a list) with self.

        If the argument is a list of Trees objects and if the topology of
        each Trees object is compatible, return a Trees object."""
        ctree_list=[]
        types=list(self.__types)
        attributes=list(self.__attributes)
        for t in range(len(tree_list)):
            if issubclass(tree_list[t].__class__, Trees):
                ctree_list.append(tree_list[t].__ctrees)
                types+=list(tree_list[t].__types)
                attributes+=list(tree_list[t].__attributes)
            else:
                ctree_list.append(tree_list[t])
        try:
            cmerged=self.__ctrees.MergeVariable(ctree_list)
        except RuntimeError, error:
            raise FormatError(error)
        merged=Trees(cmerged, types, attributes)
        self._copy_vid_conversion(merged)
        self._copy_tid_conversion(merged)
        return merged

    def MTGVertexId(self, TreeId, vid=None):
        """Return the MTG vid of a Tree vertex for a given tree"""
        if self.__tree_to_mtg_vid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if vid is None:
            return dict(self.__tree_to_mtg_vid[TreeId])
        else:
            return self.__tree_to_mtg_vid[TreeId][vid]

    def MTGComponentRoot(self, TreeId=None):
        """Return the MTG ComponentRoot corresponding to a given 
        tree identifier"""
        if self.__tree_to_mtg_tid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if TreeId is None:
            return dict(self.__tree_to_mtg_tid)
        else:
            try:
                res=self.__tree_to_mtg_tid[TreeId]
            except KeyError:
                msg="Tree number "+str(TreeId)+" of self" + \
                    "does not correspond to any component of the MTG"
                raise IndexError, msg
            else:
                return res

    def NbInt(self):
        """Return the number of variables with integer type."""
        return self.__ctrees.NbInt()

    def NbFloat(self):
        """Return the number of variables with floating type."""
        return self.__ctrees.NbFloat()

    def NbVariables(self):
        """Return the number of variables (attributes) of a tree."""
        return len(self.__types)

    def NbTrees(self):
        """Return the number of trees."""
        return self.__ctrees.NbTrees()

    def Plot(self, ViewPoint="Data", Length=None, BottomDiameter=None, 
             Color=None, DressingFile=None, Title="", variable=0):
        """Graphical output using the Geom 3D viewer for Trees 
           or Gnuplot.py for features.
        
        Usage:  Plot(ViewPoint="Data")
                Plot("FirstOccurrenceRoot", variable=0)
        Other possible values for ViewPoint: 
            FirstOccurrenceLeaves
            SojournSize
            Counting"""
        if type(ViewPoint)!=str:
            msg='bad type for argument "ViewPoint": '\
              +"type 'str' expected"
            raise TypeError, msg
        import os
        if ViewPoint.upper()=="DATA":
            # Graphical output of the Trees using the Geom 3D viewer 
            # create the MTG file
            mtgprefix="gvtmp"
            file_created=False
            while not file_created:
                try:
                    cfile=open(mtgprefix+'.mtg','r')
                except IOError:
                    # file does not exist
                    mtgfile_name= mtgprefix+".mtg"
                    file_created=True
                else:
                    cfile.close()
                    import random
                    mtgprefix+=str(random.randint(1,9))
    
            if DressingFile is None:
                # create the dressing file
                drfprefix="dftmp"
                file_created=False
                while not file_created:
                    try:
                        cfile=open(drfprefix+'.drf','r')
                    except IOError:
                        # file does not exist
                        drffile_name= drfprefix+".drf"
                        file_created=True
                    else:
                        cfile.close()
                        import random
                        drfprefix+=str(random.randint(1,9))
                
                dressing=file(drffile_name,'w+')
                dressing.write("Phyllotaxy = 103\n")
                dressing.write("NbPlantsPerLine = "+str(self.NbTrees())+"\n")
                dressing.close()
            else:
                # use the given dressing file
                drffile_name=DressingFile
    
            # create the temporary MTG file
            self.Save(mtgfile_name, False, list(self.__attributes))
    
            import openalea.aml as amlPy
            
            mode=False
            if not amlPy.getmode():
            # conversion from AML object to Python
                mode=True
                amlPy.setmode(1)
                
            M=amlPy.MTG(mtgfile_name)
            DR=amlPy.DressingData(drffile_name)
                            
            if Length is None:
                # define a default length function
                default_lengthfunc= lambda x: 20
                for var in range(self.NbVariables()):
                    lvariable_name=self.__attributes[var]
                    if ((lvariable_name.upper()).find("LEN") != -1):
                        # current variable name contains "LEN"
                        default_lengthfunc= \
                            lambda x: amlPy.Feature(x, lvariable_name)
                        break
            else:
                default_lengthfunc= Length
            # check the type and value
            def lengthfunc(x):
                res=default_lengthfunc(x)
                if (type(res)==int) or (type(res)==float):
                    if (res > 0):
                        return res
                    else:
                        return 10
                else:
                    return 10
    
            if BottomDiameter is None:
                # define a default diameter function
                default_diamfunc= lambda x: 5
                for var in range(self.NbVariables()):
                    dvariable_name=self.__attributes[var]
                    if ((dvariable_name.upper()).find("DIA") != -1):
                        # current variable name contains "DIA"
                        default_diamfunc= \
                            lambda x: amlPy.Feature(x, dvariable_name)
                        break
            else:
                default_diamfunc= BottomDiameter
            # check the type and value
            def diamfunc(x):
                res=default_diamfunc(x)
                if (type(res)==int) or (type(res)==float):
                    if (res > 0):
                        return res
                    else:
                        return 10
                else:
                    return 10
    
            if Color is None:
                # define a default color function
                default_colorfunc= lambda x: 0
                for var in range(self.NbVariables()):
                    cvariable_type=self.__types[var]
                    cvariable_name=self.__attributes[var]
                    if (cvariable_type == VariableType.STATE):
                        # current variable is a state variable
                        default_colorfunc= \
                            lambda x: amlPy.Feature(x, cvariable_name)+2
                        break
            else:
                default_colorfunc= Color
            # check the type and value
            def colorfunc(x):
                res=default_colorfunc(x)
                if (type(res)==int):
                    if (res > 0):
                        return res
                    else:
                        return 0
                else:
                    return 0
    

            vtx_list=amlPy.VtxList(Scale=1)
            pf=amlPy.PlantFrame(vtx_list, Scale=2, DressingData=DR, 
                                Length=lengthfunc, BottomDiameter=diamfunc)
            amlPy.Plot(pf, Color=colorfunc)
            
            if mode:
                amlPy.setmode(0)
            
            # remove the temporary files
            os.remove(mtgfile_name)
            if DressingFile is None:
                os.remove(drffile_name)
                
        else:
            # Graphical output of the features using Gnuplot.py
            if ViewPoint.upper()=="FIRSTOCCURRENCEROOT":
                ftype=CharacteristicType.FIRST_OCCURRENCE_ROOT
            elif ViewPoint.upper()=="FIRSTOCCURRENCELEAVES":
                ftype=CharacteristicType.FIRST_OCCURRENCE_LEAVES
            elif ViewPoint.upper()=="SOJOURNSIZE":
                ftype=CharacteristicType.SOJOURN_SIZE
            elif ViewPoint.upper()=="COUNTING":
                ftype=CharacteristicType.NB_ZONES
            else:
                msg='bad value for argument "ViewPoint": '+ViewPoint
                raise ValueError, msg
            cvariable = self._valid_cvariable(variable)
            print "Using cvariable", cvariable
            if (not self.__ctrees.IsCharacteristic(cvariable, ftype)):
                msg = "Characteristic " + ViewPoint + " not computed " + \
                    "for variable ", variable
                raise FormatError(msg)
##            ref_file_id=str(self._valid_cvariable(variable)+1)
##            # part of the filename which identifies the graph to be displayed
            file_id = str(self._valid_cvariable(variable)+1)+str(ftype+1)
##            prefix="ftmp"
##            file_created=False
##            file_list=[]
##            # find a non existing file name
##            while not file_created:
##                try:
##                    cfile=open(prefix+'11.plot','r')
##                except IOError:
##                    file_created=True
##                else:
##                    import random
##                    prefix+=str(random.randint(1,9))
            try:
                self.__ctrees.plot(Title=Title, Suffix=file_id)
                # build the list of the created files: 
##                for var in range(self.NbInt()):
##                    for char in [str(c) for c in range(5)]+[""]:
##                        filename=prefix+str(self._valid_cvariable(var)+1)+char
##                        try:
##                            tmpfile=open(filename+'.plot', 'r')
##                        except IOError:
##                            pass
##                        else:
##                            tmpfile.close()
##                            # add the .plot and .print files
##                            file_list+=[filename+extension
##                                for extension in [".plot", ".print"]]
##                            if (char=='1') or (char==''):
##                                # add the .dat file
##                                file_list+=[filename+".dat"]
            except RuntimeError, f:
##                for tmpfile in file_list:
##                   os.remove(tmpfile)
                raise FormatError, f
##            try:
##                # check if the desired file exists
##                cfile=open(prefix+file_id+'.plot','r')
##            except IOError:
##                for tmpfile in file_list:
##                   os.remove(tmpfile)
##                msg='Characteristic not computed: '+str(ViewPoint)
##                raise ValueError, msg                
            # else:
            #     nb_windows= \
            #       self.__ctrees.NbValues(self._valid_cvariable(variable))
            #     self.__plot= _PlotManager(file_list, prefix+file_id, 
            #                                          nb_windows)

    def Save(self, file_name, overwrite=False, variable_names=None):
        """Save trees into a file as a MTG."""
        if variable_names is None:
            variable_names=self.__attributes
        elif len(variable_names)!=self.NbVariables():
            msg="bad number of variable names: "+str(len(variable_names)) \
                +"; "+str(self.NbVariables())+" name(s) expected"
            raise IndexError, msg
        for var in range(len(variable_names)):
            if type(variable_names[var])!=str:
                msg="bad type for name of variable "+str(var)
                raise TypeError, msg
            elif variable_names[var].find(" ")!=-1:
                msg="name of variable "+str(var)+" should not contain "\
                     +"any space character"
                raise ValueError, msg
        if not overwrite:
            try:
                f=file(file_name, 'r')
            except IOError:
                pass # f=file(file_name, 'w+')
            else:
                msg="File "+file_name+" already exists"
                raise IOError, msg
                                
        f=file(file_name, 'w+')
        # write the MTG header
        mtgs, max_depth=self.__mtg_header(variable_names)
        f.write(mtgs)
        # print each tree
        for t in range(self.NbTrees()):
            current_tree= self.Tree(t)
            if (t > 0):
                f.write("/")
            f.write("N"+str(t)+"\n")
            f.write("\t/")
            post_tabulation="\t"*(max_depth-current_tree._ctree().Order())
            f.write(current_tree._mtg_write("\t", post_tabulation))
        f.close()
            
    def SelectVariable(self, variables, mode="Keep"):
        """Select variables.

        Select the given variable(s) and build a new Trees object.
        Argument mode must be "Keep" or "Reject" and corresponds to the
        conservation or rejection of the selected variables."""
        cvariables=[]
        variable_list=[] # list of the variables
        types=[]
        attributes=[]
        if type(variables)==int:
            # a single variable is given
            cvariables.append(self._valid_cvariable(variables)+1)
            variable_list.append(variables)
        else:
            # a list of variables is given
            for variable in variables:
                cvariables.append(self._valid_cvariable(variable)+1)
                # correspondence of variables between self and self.__ctrees
            variable_list=variables
        try:
            if not (mode.upper()=="KEEP" or mode.upper()=="REJECT"):
                msg="unknown selection mode: "+str(mode) \
                    + " - expecting 'Keep' or 'Reject'"
                raise ValueError, msg
            else:
                if mode.upper()=="KEEP":
                    keep=True
                else:
                    keep=False
                for variable in range(self.NbVariables()):
                # building the list of types and attributes
                    if (variable_list.__contains__(variable) and keep):
                        types.append(self.__types[variable])
                        attributes.append(self.__attributes[variable])
                    elif not (variable_list.__contains__(variable) or keep):
                        types.append(self.__types[variable])
                        attributes.append(self.__attributes[variable])
                cselected=self.__ctrees.SelectVariable(cvariables, keep)
        except RuntimeError, error:
            # raise stat_tool.Format_error, str(cerror)
            replaced=str(error) # self.__replacestr(str(cerror), "variable")
            raise FormatError, replaced
        selected=Trees(cselected, types, attributes)
        self._copy_vid_conversion(selected)
        self._copy_tid_conversion(selected)
        return selected

    def SelectIndividual(self, identifiers, mode="Keep"):
        """Select individuals.

        Select the given tree(s) and build a new Trees object.
        Argument mode must be "Keep" or "Reject" and corresponds to the
        conservation or rejection of the selected trees."""
        id_list=[] # list of the individuals
        if type(identifiers)==int:
            # a single identifier is given
            id_list.append(identifiers)
        else:
            # a list of identifiers is given
            id_list=list(identifiers)
        try:
            if not (mode.upper()=="KEEP" or mode.upper()=="REJECT"):
                msg="unknown selection mode: "+str(mode) \
                    + " - expecting 'Keep' or 'Reject'"
                raise ValueError, msg
            else:
                if mode.upper()=="KEEP":
                    keep=True
                else:
                    keep=False
                # build the list of types and attributes
                attributes=self.__attributes
                types=self.__types
                cselected=self.__ctrees.SelectIndividual(id_list, keep)
        except RuntimeError, error:
            raise FormatError, str(error)
        selected=Trees(cselected, types, attributes)
        self._copy_vid_conversion(selected)
        self._copy_tid_conversion(selected)
        return selected

    def SegmentationExtract(self, variable, values, mode="Keep"):
        """Return subtrees with homogeneous zones, with respect to a given
        variable and its value(s).

        The values argument can be an integer or a list of integers
        Usage: SegmentationExtract(variable, value, mode="Reject")
               SegmentationExtract(variable, values, mode="Reject")"""
        ## variable_list=[] # list of the variables
        types=[]
        attributes=[]
        cvariable=self._valid_cvariable(variable)+1
        # correspondence of variables between self and self.__ctrees
        try:
            if not (mode.upper()=="KEEP" or mode.upper()=="REJECT"):
                msg="unknown selection mode: "+str(mode) \
                    + " - expecting 'Keep' or 'Reject'"
                raise ValueError, msg
            else:
                if mode.upper()=="KEEP":
                    keep=True
                else:
                    keep=False
            for v in range(self.NbVariables()):
            # building the list of types and attributes
                if (v!=variable):
                    types.append(self.__types[v])
                    attributes.append(self.__attributes[v])
            csegmentation=self.__ctrees.SegmentationExtract(cvariable, 
                                                            values, keep)
        except RuntimeError, error:
            # raise stat_tool.Format_error, str(cerror)
            replaced=self.__replacestr(str(error), "variable") # str(error)
            raise FormatError, replaced
        segmentation=Trees(csegmentation, types, attributes)
        # self._copy_vid_conversion(selected)
        # self._copy_tid_conversion(selected)
        return segmentation

    def Shift(self, variable, param):
        """Shift all values of self corresponding to the given variable,
        from a given parameter."""
        cvariable=self._valid_cvariable(variable)+1
        if ((type(param)!=float) and (type(param)!=int)):
            msg='bad type for argument "param": '\
              +"type 'float' or 'int' expected"
            raise TypeError, msg
        # if the variable is of floating type, conversion of the parameter
        # to a float
        if ((self.Type(variable)==VariableType.REAL_VALUE) 
            and (type(param)==int)):
            eparam=param + 0.
        else:
            eparam=param
        try:
            cshifted=self.__ctrees.Shift(cvariable, eparam)
        except RuntimeError, error:
            replaced=self.__replacestr(str(error), "variable")
            raise FormatError, replaced
        shifted=Trees(cshifted, self.__types, self.__attributes)
        self._copy_vid_conversion(shifted)
        self._copy_tid_conversion(shifted)
        return shifted

    def Transcode(self, variable, new_values):
        """Transcoding of values.

        Transcode the values of a given variable.
        Argument variable refers to the transcoded variable.
        Argument new_values refers to the list of new values. """
        if ((self.Type(variable)!=VariableType.INT_VALUE) and
            (self.Type(variable)!=VariableType.STATE)):
            msg="Bad variable: " + str(variable) + ": must be of integer" +  \
                " or of state type"
            raise IndexError, msg
        cvariable=self._valid_cvariable(variable)+1
        # correspondence of variables between self and self.__ctrees
        try:
            ctranscoded=self.__ctrees.Transcode(cvariable, new_values)
        except RuntimeError, error:
            msg=str(error)
            raise FormatError, msg
        transcoded=Trees(ctranscoded, self.__types, self.__attributes)
        self._copy_vid_conversion(transcoded)
        self._copy_tid_conversion(transcoded)
        return transcoded

    def Tree(self, TreeId):
        """Return the tree corresponding to the given identifier."""
        # d=self.__default_tree_value()
        d = self.Types()
        if self.__valid_tree(TreeId):
            res = Tree(d, self.__ctrees.Tree(TreeId), self.__attributes)
            if len(self.__mtg_to_tree_vid[TreeId]) > 0:
                res._copy_vid_conversion(self.TreeVertexId(TreeId),
                                         self.MTGVertexId(TreeId))
                res._copy_mtg_tid(self.MTGComponentRoot(TreeId))
            return res

    def TreeId(self, MtgComponentRoot=None):
        """Return the tree identifier corresponding to a given MTG 
        ComponentRoot."""
        if self.__mtg_to_tree_tid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if MtgComponentRoot is None:
            return dict(self.__mtg_to_tree_tid)
        else:
            try:
                res=self.__mtg_to_tree_tid[MtgComponentRoot]
            except KeyError:
                msg="MTG Vertex "+str(MtgComponentRoot)+" is not the root "+\
                    "of any tree in self"
                raise IndexError, msg
            else:
                return res

    def TreeVertexId(self, TreeId=None, MTGVid=None):
        """Return the tree vid of a MTG vertex for a given tree, 
        or both the tree and the vid for a given MTG vertex.
        
        Usage: v = TreeVertexId(TreeId=0, MTGVid=2)
               dic = TreeVertexId(TreeId=0)
               t, v = TreeVertexId(MTGVid=2)"""
        if self.__mtg_to_tree_vid is None:
            raise Warning, "Current Trees object has not been obtained from " \
                "a MTG"
        if not(TreeId is None):
            if MTGVid is None:
                return dict(self.__mtg_to_tree_vid[TreeId])
            else:
                return self.__mtg_to_tree_vid[TreeId][MTGVid]
        else:
            if MTGVid is None:
                raise ValueError, "TreeVertexId requires at least one argument"
            else:
                # find MTGVid in all the dictionaries
                correct_vid=False
                for t in range(self.NbTrees()):
                    try:
                        tree_vid=self.__mtg_to_tree_vid[t][MTGVid]
                    except KeyError:
                        pass
                    else:
                        correct_vid=True
                        break
                if not(correct_vid):
                    msg = "MTG Vertex "+str(MTGVid)+" not found"
                    raise ValueError, msg
                else:
                    return t, tree_vid

    def Type(self, index):
        """Return the type of one given variable (attribute).

        Argument index, which must be in the range [0, NbVariables()],
        refers to the concerned variable."""
        return self.__types[index]
        # return self.__ctrees.Type(self.__tmap[index])

    def Types(self):
        """Return the type of each tree variable (attribute)."""
        return list(self.__types)

    def _ctrees(self):
        return self.__ctrees

    def _ctrees_display(self):
        return self.__ctrees.Display(False)
    
    def _copy_vid_conversion(self, dest):
        # copy the dictionnaries corresponding to the tree -> MTG
        # and MTG -> tree vid conversion
            if self.__mtg_to_tree_vid is None:
                dest.__mtg_to_tree_vid=None
            else:
                dest.__mtg_to_tree_vid=list(self.__mtg_to_tree_vid)
            if self.__tree_to_mtg_vid is None:
                dest.__tree_to_mtg_vid=None
            else:
                dest.__tree_to_mtg_vid=list(self.__tree_to_mtg_vid)

    def _copy_tid_conversion(self, dest):
        # copy the dictionnaries corresponding to the tree -> MTG
        # and MTG -> tree id conversion
            if self.__mtg_to_tree_tid is None:
                dest.__mtg_to_tree_tid=None
            else:
                dest.__mtg_to_tree_tid=dict(self.__mtg_to_tree_tid)
            if self.__tree_to_mtg_tid is None:
                dest.__tree_to_mtg_tid=None
            else:
                dest.__tree_to_mtg_tid=dict(self.__tree_to_mtg_tid)

    def _max(self, variable=None):
        """Return the maximal values of one or all variables."""
        i=self.__ctrees.Max()
        for v in range(len(i)):
            if self.__types[v]==VariableType.REAL_VALUE:
                i[v]+=0.
        res=TreeValue(i)
        if variable is None:
            return res
        else:
            self._valid_cvariable(variable)
            return res[variable]

    def _min(self, variable=None):
        """Return the minimal values of one or all variables."""
        i=self.__ctrees.Min()
        for v in range(len(i)):
            if self.__types[v]==VariableType.REAL_VALUE:
                i[v]+=0.
        res=TreeValue(i)
        if variable is None:
            return res
        else:
            self._valid_cvariable(variable)
            return res[variable]
    
    def _set_int_value_type(self, variable):
        # turns a state variable into an integer-valued variable
        if (variable < 0 or variable >= self.NbVariables()):
            raise IndexError, "variable index out of range: "+str(variable)
        elif (self.__types[variable]==VariableType.STATE):
            self.__types[variable]=VariableType.INT_VALUE
        else:
            msg="bad type for variable " + str(variable) + ": " + \
                str(self.__types[variable])
            raise TypeError, msg
    
    def _valid_cvariable(self, variable):
        # check the validity of the variable index and return the corresponding 
        # index for the C++ programs
        if (type(variable)!=int):
            raise TypeError, "bad variable index type: "+str(type(variable))
        elif (variable < 0 or variable >= self.NbVariables()):
            raise IndexError, "variable index out of range: "+str(variable)
        else:
            return self.__tmap[variable]

    def __default_tree_value(self):
        l=[]
        for var in range(self.NbVariables()):
            if self.__types[var] == VariableType.REAL_VALUE:
                l.append(0.)
            else:
                l.append(0)
        default_value=TreeValue(l)
        return default_value
    
    def __mtg_header(self, variable_names):
        # return the MTG header and the maximal tree depth
        max_depth=0
        # computation of the maximal depth in any tree
        for t in range(self.NbTrees()):
            max_depth=max(max_depth, self.Tree(t)._ctree().Order())

        mtgs='#one-scaled MTG generated by module "trees"\n' 
        mtgs+="CODE:\tFORM-A\n\nCLASSES:\n"
        mtgs+="SYMBOL\tSCALE\tDECOMPOSITION\tINDEXATION\tDEFINITION\n"
        mtgs+="$\t0\tFREE      \tFREE      \tIMPLICIT\n"
        mtgs+="N\t1\tCONNECTED \tFREE      \tIMPLICIT\n"
        mtgs+="V\t2\tNONE      \tFREE      \tEXPLICIT\n\n"
        mtgs+="DESCRIPTION:\n"
        mtgs+="LEFT\tRIGHT\tRELTYPE\tMAX\n"
        mtgs+="V\tV\t+\t?\n"
        mtgs+="V\tV\t<\t?\n\n"
        mtgs+="FEATURES:\n"
        mtgs+="NAME\tTYPE\n\n"
        for var in range(self.NbVariables()):
            mtgs+=variable_names[var]+"\t"
            if self.Type(var)==VariableType.REAL_VALUE:
                mtgs+="REAL"
            else:
                mtgs+="INT"
            mtgs+="\n"
        mtgs+="\nMTG:\nENTITY-CODE"
        mtgs+="\t"*(max_depth+1)
        for var in range(self.NbVariables()):
            mtgs+=variable_names[var]+"\t"
        mtgs+="\n/"
        return [mtgs, max_depth]
    
    def __name_to_variable(self, name):
        # convert a variable name into a variable number
        if type(name)!=str:
            msg="bad type for variable name: "+str(type(name)) \
                + " - excpeting type 'str'"
            raise TypeError, msg
        variable=0
        found=False
        while ((variable < self.NbVariables()) and (not found)):
            if name.upper() == self.__attributes[variable].upper():
                found=True
            else:
                variable+=1
        if not found:
            msg="bad variable name: "+str(name)
            raise ValueError, msg
        else:
            return variable

    def __replacestr(self, message, sub, upcase=False):
        # applies the variable mapping defined by self.__tmap
        # to string messages
        if upcase:
            sub=sub.upper()
        index=0
        while index < len(message):
            i=message.find(sub, index)
            if (i==-1) or (i+len(sub)+1 >= len(message)):
                # subchain has not been found or this is the last word
                return message
            else:
                index=i+1
                sval=message[i+len(sub)+1]
                # caracter following sub
                try:
                    val=int(sval)
                except ValueError:
                    pass
                else:
                    message=message[0:i+len(sub)+1]+str(self.__tmapi[val-1]) \
                            +message[i+len(sub)+2:len(message)]
        return message
    
    def __valid_tree(self, TreeId):
        if type(TreeId)!=int:
            raise TypeError, "bad tree index type: "+str(type(TreeId))
        elif (TreeId < 0) or (TreeId >= self.NbTrees()):
            raise IndexError, "tree index out of range: "+str(TreeId)
        else:
            return True
        
    def __str__(self):
        """Display the trees."""
        
        # replaced=self.__replacestr(self.__ctrees.Display(False), 
        #                            "variable", True)
        # return replaced
        classstr=str(self.__class__)
        res=classstr+": "+str(self._ctrees())
        return res


if __name__ == '__main__':
    pass # add a call to run your script here


if __name__ == '__main__':
    pass # add a call to run your script here


if __name__ == '__main__':
    pass # add a call to run your script here
