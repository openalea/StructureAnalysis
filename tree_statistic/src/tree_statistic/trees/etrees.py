"""Trees with integral and floating attributes with enhanced methods.
"""


import trees, ctree, int_fl_containers
import openalea.stat_tool

stat_tool=openalea.stat_tool
I_DEFAULT_TREE_SIZE=trees.I_DEFAULT_TREE_SIZE
I_DEFAULT_TREE_DEPTH=trees.I_DEFAULT_TREE_DEPTH
# VariableType=trees.VariableType

TreeValue=trees.TreeValue
    
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
        Examples: T=Tree([0], 0, 0)
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
                if t==stat_tool.VariableType.REAL_VALUE:
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
        referred to by its identifier (vid)."""
        if vid is None:
            return trees.Tree.Depth(self)
        else:
            self.__valid_vid(vid)
            return self.__ctree.Depth(vid)
    
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
       
    def NbChildren(self, vid):
        """Return the number of children of a given vertex,
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.NbChildren(vid)

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
        elif issubclass(arg.__class__, stat_tool.Distribution):
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
       
    def NbChildren(self, vid):
        """Return the number of children of a given vertex,
        referred to by its identifier (vid)."""
        self.__valid_vid(vid)
        return self.__ctree.NbChildren(vid)

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
