#
#
#

"""

"""

def signat_vert_comp(g):
    """
    """
    signat_vertices={}
    for vid in g.vertices():
        s=[]
        for vid2 in g.vertices():
            if vid2 in g.out_neighbors(vid):
                for eid in g.edges():
                    if g._edges[eid]==(vid,vid2):
                        s.append(g._edge_property['weight'][eid])
            else:
                s.append(0)
        signat_vertices[vid]=s                             
    return signat_vertices

def signat_edg_comp(g, signat_vertices):
    '''    
    '''
    signat_edges={}
    for eid in g.edges():
        vid1=g._edges[eid][0]
        vid2=g._edges[eid][1]
        s=signat_vertices[vid1][:]       
        s.append(s[vid2-1])
        s[vid2-1]=0
        signat_edges[eid]=s
    return signat_edges

def quasi_isomorphic_edges(eid1, eid2, signat_edges):
    '''
    '''
    if signat_edges[eid1]==signat_edges[eid2]:
        return True
    else:
        return False

def Rre_compt(g, paths):
    '''Computation of the reduction graph with return edges.

    Each path will be transformed into a loop...
    
    :Parameters:
      - `g` (:class:`RootedGraph`) - The graph to be reduced
      - `paths` (list of list of vertices) - List of Quasi isomorph vertices.

    :Returns:
      - The initial graph transformed into a cyclic graph.

    :Examples:

    >>> Rre_compt(g,[[1,2,3], [4,5,6]])

    .. note:: The return graph is a modification of the initial one.
    .. seealso:: :class:`RootedGraph`
    '''

    return g,

def C_signat_edg(edge_signature):
    '''Computation of the reduction graph with return edges.

    Each path will be transformed into a loop...
    
    :Parameters:
      - `edge_signature` (dict{eid:[int]}) - Signature of each edge

    :Returns:
      - Set of edge signature values.

    :Examples:

    >>> Rre_compt(g,[[1,2,3], [4,5,6]])

    .. note:: The return graph is a modification of the initial one.
    .. seealso:: :class:`RootedGraph`    
    '''
    return set(edge_signature.values())

def overlapped_QPP_sets(in1):
    '''    
    '''

    return [out1]

def QIE_sets_compt(in1, in2):
    '''    
    '''

    return [out1]

def QPPs_compt(in1):
    '''    
    '''
  
    return [out1]

def selected_QPP(in1, in2):
    '''    
    '''

    return [out1]

def rooted_graph():
    '''    
    '''

    return [out1]

def rooted_graph_representation(in1):
    '''    
    '''
    
    return []

