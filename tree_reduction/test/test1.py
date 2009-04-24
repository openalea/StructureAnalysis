# -*- coding: cp1252 -*-
from vplants.tree_reduction.graph import *
   
def test1():
    edges=[(1,2),(2,3),(3,4),(3,5),(3,6)]
    root = 1
    g = from_edges(root, edges)
    assert g.nb_vertices() == 6
    assert g.nb_edges() == 5
    return g
    
