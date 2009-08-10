# -*- coding: utf-8 -*-

'''
This class is a temporary class.
It will be replaced by a rooted graph class in openalea.container.
'''

from openalea.container.property_graph import PropertyGraph

class RootedGraph(PropertyGraph):
    def __init__(self, root=0, graph=None):
        PropertyGraph.__init__(self,graph)
        self._root = root
        self.add_edge_property('weight')        
        self.add_edge_property('loop')
        
    def get_root(self):
        return self._root
    def set_root(self, root):
        self._root = root
    root = property(get_root, set_root)

def from_edges(root, edges, edges_weight):
    g = RootedGraph(root=root)
    for e in edges:
        s,t=e
        if s not in g:
            g.add_vertex(s)
        if t not in g:
            g.add_vertex(t)
        g.add_edge(s,t)
    weight = g.edge_property('weight')
    for eid in g.edges():
        weight[eid]=edges_weight[eid]  
    update_loop(g)
    return g

def update_loop(g):
    loop = g.edge_property('loop')
    for eid in g.edges():
        loop[eid]=-1
        
def topological_sort(g, vtx_id):
    ''' 
    Traverse the directed graph `g` in a prefix way.
    (root then children)

    This is a non recursive implementation.
    '''
    yield vtx_id
    for vid in g.out_neighbors(vtx_id):
        for node in topological_sort(tree, vid):
            yield node

    
