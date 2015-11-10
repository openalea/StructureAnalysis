# -*- coding: utf-8 -*-
"""Trees with integral and floating attributes with enhanced methods.
"""


import trees, ctree
import openalea.tree_statistic.int_fl_containers as int_fl_containers
import openalea.tree_statistic._errors as _errors
import openalea.stat_tool as stat_tool
import openalea.stat_tool.error as check_error

VariableType = stat_tool.VariableTypeBis

I_DEFAULT_TREE_SIZE = trees.I_DEFAULT_TREE_SIZE
I_DEFAULT_TREE_DEPTH = trees.I_DEFAULT_TREE_DEPTH

StatTreeError = _errors.StatTreeError
VariableTypeDict = VariableType.values

TreeValue = trees.TreeValue


class Tree(trees.Tree):
    """An implementation of trees with enhanced functionalities."""
    
    def __init__(self, arg, arg2=None, root=None):
        """Initialize a Tree from a list of values or a tree object.
        
        Initialized Tree object can be defined either from: 
        - a list of integral and floating values or a TreeValue object,
          the tree size and the root identifier;
        - a list of integral and floating values or a TreeValue object
          and a ctree.CTree object;
        - a list of integral and floating values or a TreeValue object
          and a TreeStructure object;
        - a Tree object
        Examples: T = Tree([0], 0, 0)
                  U=Tree([0, 1.], TreeStructure(T))""" 
        if issubclass(arg.__class__, trees.Tree):
            #arg is supposed to be a Tree...            
            trees.Tree.__init__(self, arg)
            # self.__types=arg.Types()
        else:
            #... or a list of values
            default_value=TreeValue(arg)
            self.__types=list(default_value.Types())
            if issubclass(arg2.__class__, ctree.CTree):
                # arg2 is supposed to be a CTree...
                self.__ctree=ctree.CTree(arg2)
            elif issubclass(arg2.__class__, trees.TreeStructure):
                # ... or a tree structure
                trees.Tree.__init__(self, arg, arg2)
            else:
                #... or the number of vertices
                if arg2 is None:
                    size=0
                else:
                    size=arg2
                if root is None:
                    troot=0
                else:
                    troot=root
                self.__ctree=ctree.CTree(default_value.NbInt(), 
                                          default_value.NbFloat(), size, troot)
            # complete the fields related to MTGs
            self.__mtg_to_tree_vid=None
            self.__tree_to_mtg_vid=None
            self.__mtg_tid=None
            attributes=[]
            for var in range(self.NbVariables()):
                attributes.append("Variable"+str(var))
            self.__attributes=list(attributes)

    def AddVertex(self, value=None):
        """Add a vertex with given value and return its identifier.
        
        Argument value can be a TreeValue, a list that can be converted
        to a TreeValue, or it can be omitted. In the latter case, null values
        will be used."""
        
        if value is None:
            lvalue=[]
            for t in self.__types:
                if t==VariableType.REAL_VALUE:
                    lvalue.append(0.)
                else:
                    lvalue.append(0)
            tvalue=TreeValue(lvalue)
        else:
            tvalue=TreeValue(value)
        return self.__ctree.AddVertex(tvalue)

    def AddEdge(self, parentvid, childvid, edgetype='+'):
       """Add an edge of given type between two vertices referred to by their 
       vertex identifiers (vids)."""
       self.__valid_vid(parentvid)
       self.__valid_vid(childvid)
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
       self.__ctree.AddEdge(parentvid, childvid, btype)

    def Children(self, vid):
        """Get an iterator on the children of a given vertex
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.Children(vid)

    def Depth(self, vid=None):
        """Return tree depth or the depth of a given vertex,
        referred to by its identifier (vid).

        The depth of a vertex is the number of vertices along the path
        between this vertex and the root vertex
        (including both vertices, which may be confounded).
        The depth of a tree is the maximal depth of its vertices."""
        if vid is None:
            return trees.Tree.Depth(self) + 1
        else:
            self.__valid_vid(vid)
            return self.__ctree.Depth(vid) + 1

    def Display(self, vids=True, attributes=True, key=None, mtg_vids=False):
        """Display the subtree rooted at argument key, with or without 
        the vertex identifiers (vids) and the attributes.
        
        The attributes are displayed (after the vids) if and only if 
        boolean argument attributes is True. If attributes is False, only the
        vids are displayed, else the vids are displayed if and only if 
        vids is True."""
        if key is None:
            key=self.Root()
        if not (vids or attributes):
            vids=True
        print self._display(key, vids, attributes, mtg_vids)
    
    def IsRoot(self, vid):
        """Return True if and only if argument is the tree root 
        vertex identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.IsRoot(vid)
       
    def NbChildren(self, vid, EdgeType=None):
        """Return the number of children of a given vertex,
        referred to by its identifier (vid). Filter according to EdgeType"""
        self.__valid_vid(vid)
        if (EdgeType is None):
            return self.__ctree.NbChildren(vid)
        else:
            c = list(self.Children())
            return len([v for v in c if self.EdgeType(vid, v) == EdgeType])

    def Parent(self, vid):
        """Return the parent of the given vertex,
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.Parent(vid)

    # iterators 
    def Breadthorder(self):
        """Get an iterator on the tree vertices, using a breadth-first 
        tree traversing."""
        return self.__ctree.Breadthorder()
        
    def Inorder(self):
        """Get an iterator on the tree vertices, using an inorder traversing."""
        return self.__ctree.Inorder()

    def LeavesFirstorder(self):
        """Get an iterator on the tree vertices. The leaf nodes come first, 
        then their parents, etc."""
        return self.__ctree.LeavesFirstorder()
    
    def Postorder(self):
        """Get an iterator on the tree vertices, using a postorder traversing."""
        return self.__ctree.Postorder()
    
    def Preorder(self):
        """Get an iterator on the tree vertices, using a preorder traversing."""
        return self.__ctree.Preorder()
    
    def __iter__(self):
        """Get an iterator on the tree vertices."""
        return self.__ctree.Vertices()  

class TreeStructure(trees.TreeStructure):
    """An implementation of tree structures with enhanced functionalities."""
    
    def __init__(self, arg=None, arg2=I_DEFAULT_TREE_SIZE, 
                 arg3=None):
        """Initialize a TreeStructure from a tree object, a distribution 
        or the number or vertices.
    
        Initialize a TreeStructure object from either: 
        - the number of vertices and the root identifier;
        - a Tree or a TreeStructure object;
        - a distribution and values for the maximal size and depth.""" 
        if arg is None:
            trees.TreeStructure.__init__(self)
        elif issubclass(arg.__class__, Tree):
            # arg is supposed to be a tree... 
            trees.TreeStructure.__init__(self, arg)
        elif issubclass(arg.__class__, TreeStructure):
            #... or a tree structure
            trees.TreeStructure.__init__(self, arg)
        elif issubclass(arg.__class__, stat_tool._stat_tool._Distribution):
            #... or a Distribution
            trees.TreeStructure.__init__(self, arg, arg2, arg3)
        else:
            #... or the number of vertices
            self.__tree=ctree.TreeStructure(arg, arg2)
            
    def AddVertex(self):
        """Add a vertex and return its identifier."""        
        return self.__tree.AddVertex()

    def AddEdge(self, parentvid, childvid, edgetype='+'):
       """Add an edge between two vertices referred to by their 
       vertex identifiers (vids)."""
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
       self.__tree.AddEdge(parentvid, childvid, btype)

    def Children(self, vid):
        """Get an iterator on the children of a given vertex
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__tree.Children(vid)
    
    def Depth(self, vid=None):
        """Return tree depth or the depth of a given vertex,
        referred to by its identifier (vid)."""
        if vid is None:
            return trees.TreeStructure.Depth(self)
        else:
            self.__valid_vid(vid)
            return self.__ctree.Depth(vid)
    
    def Display(self, vids=True, attributes=True, key=None):
        """Display the subtree rooted at argument key, with or without 
        the vertex identifiers (vids) and the attributes.
        
        The attributes are displayed (after the vids) if and only if 
        boolean argument attributes is True. If attributes is False, only the
        vids are displayed, else the vids are displayed if and only if 
        vids is True."""
        if key is None:
            key=self.Root()
        if not (vids or attributes):
            vids=True
        return self.__display(key, vids, attributes)
    
    def IsRoot(self, vid):
        """Return True if and only if argument is the tree root 
        vertex identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.IsRoot(vid)
       
    def NbChildren(self, vid, EdgeType=None):
        """Return the number of children of a given vertex,
        referred to by its identifier (vid). Filter according to EdgeType"""
        self.__valid_vid(vid)
        if (EdgeType is None):
            return self.__ctree.NbChildren(vid)
        else:
            c = list(self.Children())
            return len([v for v in c if self.EdgeType(vid, v) == EdgeType])

    def Parent(self, vid):
        """Return the parent of the given vertex,
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.Parent(vid)

    # iterators 
    
    def Breadthorder(self):
        """Get an iterator on the tree vertices, using a breadth-first 
        tree traversing."""
        return self.__ctree.Breadthorder()
        
    def Inorder(self):
        """Get an iterator on the tree vertices, using an inorder traversing."""
        return self.__ctree.Inorder()
    
    def Postorder(self):
        """Get an iterator on the tree vertices, using a postorder traversing."""
        return self.__ctree.Postorder()
    
    def Preorder(self):
        """Get an iterator on the tree vertices, using a preorder traversing."""
        return self.__ctree.Preorder()

    def _SetMTGVidDictionary(self, VidDict, TreeId=None, ValidityCheck=False):
        """Set the dictionaries corresponding to the tree -> MTG
            and MTG -> tree vertex identifiers correspondences.

        :Usage:

            _SetMTGVidDictionary(self, VidDict, TreeId=None, ValidityCheck=False)

        :Parameters:

          `VidDict` (list or dict) - Dictionary or list of dictionaries with the vertices in self
            as keys and the vids of a MTG as values
          `TreeId` (int) - Identifier of the tree whose MTG ids must be set
            (all trees in self if None)
          `ValidityCheck` (bool) - Check whether the dictionary values correspond to valid 
            tree vertex identifiers

        :Remarks:

            If TreeId is None, VidDict must be a list of dictionaries with length self.NbTrees().
            Otherwise, VidDict must be a single dictionary.
            The keys of the dictionary(ies) are the vertex identifiers of the MTG,
            and values are the corresponding vertex identifiers in associated trees.
        """
        msg = "Correspondence between MTG and tree vertex identifiers "
        msg += "was previously defined already. This will be overwritten."
        if ((TreeId is None) and not(self.__mtg_to_tree_vid is None)
            and (len(self.__mtg_to_tree_vid) > 0)):
            warnings.warn(msg, Warning)
        if self.__mtg_to_tree_vid is None:
            self.__mtg_to_tree_vid = []
            for t in range(self.NbTrees()):
                self.__mtg_to_tree_vid.append({})
        if self.__tree_to_mtg_vid is None:
            self.__tree_to_mtg_vid = []
            for t in range(self.NbTrees()):
                self.__tree_to_mtg_vid.append({})
        if not(TreeId is None):
            check = self._valid_tree(Treeid)
            CpVidDict = dict(VidDict)
            VidDict = []
            for t in range(self.NbTrees()):
                if (TreeId == t):
                    VidDict.append(CpVidDict)
                else:
                    VidDict.append({})
            if (len(self.__mtg_to_tree_vid[t]) > 0):
                warnings.warn(msg, Warning)
        elif (len(VidDict) != self.NbTrees()):
            msg = "Bad number of dictionaries: " + str(len(VidDict))
            msg += " - should be " + str(self.NbTrees())
            raise ValueError, msg
        if (ValidityCheck):
            for t in range(self.NbTrees()):
                if ((TreeId is None) or (TreeId == t)):
                    for v in VidDict[t].values():
                        check = self._valid_vid(t, v)
        for t in range(self.NbTrees()):
            if ((TreeId is None) or (TreeId == t)):
                self.__mtg_to_tree_vid[t] = dict(VidDict[t])
                self.__tree_to_mtg_vid[t] = {}
                for k in VidDict[t].keys():
                    v = VidDict[t][k]
                    if self.__tree_to_mtg_vid[t].has_key(v):
                        msg = "Tree vertex " + str(v)
                        msg += " already present in dictionary for "
                        msg += "tree " + str(t)
                        raise ValueError, msg
                    else:
                        self.__tree_to_mtg_vid[t][v] = k

    def __valid_edge(self, parent, child):
        self.__valid_vid(parent)
        self.__valid_vid(child)
        valid=self.__tree.IsEdge(parent, child)
        if not(valid):
            msg="pair ("+str(parent)+", "+str(child)+") is not a valid edge"
            raise IndexError, msg
    
    def __iter__(self):
        """Get an iterator on the tree vertices."""
        return self.Vertices()

class Trees(trees.Trees):
    """An implementation of trees with enhanced functionalities."""

    def __init__(self, arg, arg2=None,
                 attribute_names=None, attribute_def=None, scale=None):
        """Initialize a Trees object from trees.

        Initialize a Trees object from either:
        - a list of Tree objects;
        - a Trees object.
        - a MTG file, a filter on the vertices, a list of attribute names,
          a list of attribute functions and the considered scale."""
        super(Trees, self).__init__(arg, arg2, attribute_names,
                                    attribute_def, scale)

    def _SetMTGVidDictionary(self, VidDict, TreeId=None, ValidityCheck=False):
        """Set the dictionaries corresponding to the tree -> MTG
            and MTG -> tree vertex identifiers correspondences.

        :Usage:

            _SetMTGVidDictionary(self, VidDict, TreeId=None, ValidityCheck=False)

        :Parameters:

          `VidDict` (list or dict) - Dictionary or list of dictionaries with the vertices in self
            as keys and the vids of a MTG as values
          `TreeId` (int) - Identifier of the tree whose MTG ids must be set
            (all trees in self if None)
          `ValidityCheck` (bool) - Check whether the dictionary values correspond to valid
            tree vertex identifiers

        :Remarks:

            If TreeId is None, VidDict must be a list of dictionaries with length self.NbTrees().
            Otherwise, VidDict must be a single dictionary.
            The keys of the dictionary(ies) are the vertex identifiers of the MTG,
            and values are the corresponding vertex identifiers in associated trees.
            
        :Examples:

        .. doctest::
            :options: +SKIP

            >>> _SetMTGVidDictionary(self, VidDict, TreeId=None, ValidityCheck=False)

        .. seealso::
            :func:`~openalea.tree_statistic.trees.Trees.MTGVertexId`,
            :func:`~openalea.tree_statistic.trees.Trees.TreeVertexId`.
        """
        msg = "Correspondence between MTG and tree vertex identifiers "
        msg += "was previously defined already. This correspondence "
        msg += "will be overwritten."
        import warnings
        replace = True
        if ((TreeId is None) and not(self.__mtg_to_tree_vid is None)
            and (len(self.__mtg_to_tree_vid) > 0)):
            warnings.warn(msg, Warning)
        if self.__mtg_to_tree_vid is None:
            self.__mtg_to_tree_vid = []
            for t in range(self.NbTrees()):
                self.__mtg_to_tree_vid.append({})
        if self.__tree_to_mtg_vid is None:
            self.__tree_to_mtg_vid = []
            for t in range(self.NbTrees()):
                self.__tree_to_mtg_vid.append({})
        if self.__tree_to_mtg_tid is None:
            self.__tree_to_mtg_tid = {}
        if self.__mtg_to_tree_tid is None:
            self.__mtg_to_tree_tid = {}
        if not(TreeId is None):
            check = self._valid_tree(Treeid)
            CpVidDict = dict(VidDict)
            VidDict = []
            for t in range(self.NbTrees()):
                if (TreeId == t):
                    VidDict.append(CpVidDict)
                else:
                    VidDict.append({})
            if (len(self.__mtg_to_tree_vid[t]) > 0):
                warnings.warn(msg, Warning)
        elif (len(VidDict) != self.NbTrees()):
            if len(VidDict) == 0:
                # erase dictionaries
                replace = False
                self.__mtg_to_tree_vid = None
                self.__tree_to_mtg_vid = None
                self.__tree_to_mtg_tid = None
                self.__mtg_to_tree_tid = None
            else:
                msg = "Bad number of dictionaries: " + str(len(VidDict))
                msg += " - should be " + str(self.NbTrees())
                raise ValueError, msg
        if (ValidityCheck):
            for t in range(self.NbTrees()):
                if ((TreeId is None) or (TreeId == t)):
                    for v in VidDict[t].values():
                        check = self._valid_vid(t, v)
                    for k in VidDict[t].keys():
                        check_error.CheckType([k], [int])
        for t in range(self.NbTrees()):
            # copy dictionary MTG->Tree
            if (((TreeId is None) or (TreeId == t)) and replace):
                self.__mtg_to_tree_vid[t] = dict(VidDict[t])
                # build dictionary Tree->MTG
                self.__tree_to_mtg_vid[t] = {}
                for k in VidDict[t].keys():
                    v = VidDict[t][k]
                    if self.__tree_to_mtg_vid[t].has_key(v):
                        msg = "Tree vertex " + str(v)
                        msg += " already present in dictionary for "
                        msg += "tree " + str(t)
                        raise ValueError, msg
                    else:
                        self.__tree_to_mtg_vid[t][v] = k
                # update dictionaries MTGComponentRoot <--> Tree Roots
                tr = self._ctrees().Tree(t).Root() # tree root
                try:
                    v = self.__tree_to_mtg_vid[t][tr]  # MTGComponentRoot
                except KeyError, error:
                    if (ValidityCheck):
                        raise KeyError, error
                    else:
                        v = sorted(VidDict[t].keys())[0]
                self.__mtg_to_tree_tid[v] = t
                self.__tree_to_mtg_tid[t] = v

def Success1stCheck(T):
    """
    Check whether < is always the first child in T

    :Parameters:
        `T` (trees.Trees) - Tree to be tested.

    :Returns:
        True if and only if < is always the first child in T
    """
    import openalea.tree_statistic.trees.etrees as etrees
    ET = etrees.Tree(T)
    for v in ET.Breadthorder():
        # get children of every vertex
        c = list(ET.Children(v))
        has_successor = False
        # types of children
        ctypes = [ET.EdgeType(v, d) for d in c]
        if (('<' in ctypes) and (ctypes[0] != '<')):
            return False
    return True

def Success1stBuild(T, Reverse=False):
    """
    Build a copy of T where < is always the first child.
    Reverse or not the other children.

    :Parameters:
        `T` (trees.Trees) - Tree to rearrange.
        `Reverse` (bool) - True iif the order of + children
            should be inversed

    :Returns:
        A pair composed by a :ref:`openalea.tree_statistic.trees.etrees` instance
        and by a dictionary d with correspondance between the vertices in T is returned
        (use as: d[T::Vertex] = Success1stBuild(T)::Vertex)

    :Examples:

    .. doctest::

        >>> T = Success1stBuild(T, Reverse=False)

    """
    import openalea.tree_statistic.trees.etrees as etrees
    ET = etrees.Tree(T)
    copy_vid_dic = True
    nvid_to_ovid = {} # conversion from new to old vids
    ovid_to_nvid = {} # conversion from new to old vids
    # conversion from Tree to MTG if any
    try:
        mtg_to_tree_vid_old = T.TreeVertex()
        tree_to_mtg_vid_old = T.MTGVertex()
        mtg_tid = T._mtg_tid() # root of mtg
    except Warning:
        copy_vid_dic = False
    # censored += [{}]
    mtg_to_tree_vid_new = {}
    tree_to_mtg_vid_new = {}
    # build a new etrees.Tree
    R = etrees.Tree(ET.Get(0), 0, 0)
    children_list = [] # list of list of children
    type_list = [] # list of list of edge types
    for v in ET.Breadthorder():
        # get children of every vertex
        c = []
        t = []
        for d in ET.Children(v):
            if ET.EdgeType(v, d) == "<":
                c = [d] + c
                t = ["<"] + t
            elif Reverse:
                c = [d] + c
                t = ["+"] + t
            else:
                c += [d]
                t += ["+"]
        children_list.append(c)
        type_list.append(t)
    lindex = 0 # index in children_ and type_lists
    # add root vertex
    v = ET.Root()
    val = ET.Get(v)
    vid = R.AddVertex(val)
    tree_to_mtg_vid_new[vid] = ET.MTGVertex(v)
    mtg_to_tree_vid_new[tree_to_mtg_vid_new[vid]] = vid
    ovid_to_nvid[v] = vid
    nvid_to_ovid[vid] = v
    for v in ET.Breadthorder():
        # add the children of v and their attributes
        llindex = 0 # index in children_ and type_[lindex] lists
        for d in children_list[lindex]:
            val = ET.Get(d)
            vid = R.AddVertex(val)
            tree_to_mtg_vid_new[vid] = ET.MTGVertex(d)
            mtg_to_tree_vid_new[tree_to_mtg_vid_new[vid]] = vid
            ovid_to_nvid[d] = vid
            nvid_to_ovid[vid] = d
            R.AddEdge(ovid_to_nvid[v], vid, type_list[lindex][llindex])
            llindex += 1
        lindex += 1
    if copy_vid_dic:
        R._copy_vid_conversion(mtg_to_tree_vid_new, tree_to_mtg_vid_new)
        R._copy_mtg_tid(mtg_tid) # copy root of mtg
    return R, ovid_to_nvid
    return R, ovid_to_nvid
