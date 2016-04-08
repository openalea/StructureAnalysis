#---------------------------------------
#Transformation of a MTG from newmtg en Type TreeGraph
#---------------------------------------


from openalea.mtg.algo import descendants
from openalea.tree_matching import TreeGraph


def mtg2treegraph(mtg,scale = 1, root = None):
    """ Return conversion of mtg into a TreeGraph structure of tree_matching.
    And corespondence map betwee, id of mtg and id of tree graph """
    tree = TreeGraph()
    
    if root is None: 
        # take first root of mtg at the given scale
        root = mtg.components_at_scale_iter(mtg.root,scale=scale).next()
    
    if mtg.scale(root) != scale:
        root = mtg.components_at_scale_iter(root,scale=scale).next() 
    
    tree.addNode(0,-1)
    idmap = { root : 0 }
    id = 1
    
    for n in descendants(mtg,root):
        if n == root: continue
        fid = mtg.parent(n)
        if fid : fid = idmap[fid]
        else: fid = -1
        tree.addNode(id,fid)
        idmap[n] = id
        id += 1
    
    return tree, idmap





