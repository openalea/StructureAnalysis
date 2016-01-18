from openalea.plantgl.scenegraph import BoundingBox, Point3Array, Polyline
from openalea.mtg.mtg import colored_tree
from util import readfile, writefile
import os, os.path

def get_shared_data(x) : return x

def compute_global_bbx(bbxprop):
    firstbbx = bbxprop.itervalues().next()
    globalbbox = BoundingBox(firstbbx.lowerLeftCorner,firstbbx.upperRightCorner)
    for bbx in bbxprop.itervalues(): globalbbox.extend(bbx)
    return globalbbox

def comparable_mtg(mtg):
    """ Create a macroscopic scale (2) with nodes corresponding to non branching parts. 
        macro nodes has associated skeleton, boundingbox and size.
        Scale 3 correspond to normal nodes. """
    s = mtg.max_scale()
    root = mtg.roots_iter(scale=s).next()  # starting node ID of scale=s, eg, in scale 3 --> starting node ID is 3

    colors = {}
    colors[1] = [root]
    colors[2] = [vid for vid in mtg.vertices_iter(scale=s) if vid == root or mtg.nb_children(mtg.parent(vid)) > 1]
    colors[3] = mtg.vertices(scale=s)

    new_mtg, mapping_ids = colored_tree(mtg, colors)
    
    positions = new_mtg.property('position')
    radii = new_mtg.property('radius')
        
    skeletons = {}
    bbox = {}
    sizes = {}
    for vid in new_mtg.vertices_iter(scale=2):
        component_iter = new_mtg.components_iter(vid)
        fcomponent = component_iter.next()
        pfcomponent = new_mtg.parent(fcomponent)
        points = []
        lradii = []
        # append position of the parent of the first element to have first segment
        if not pfcomponent is None:
            points.append(positions[pfcomponent])
            lradii.append(radii[pfcomponent])
        points.append(positions[fcomponent])
        lradii.append(radii[fcomponent])
        othercomponents = list(component_iter)
        points += [positions[cid] for cid in othercomponents]
        lradii += [radii[cid] for cid in othercomponents]
        assert len(points) > 1
        skel = Point3Array(points)
        skeletons[vid] = skel
        polskel = Polyline(skel)
        bbox[vid] = BoundingBox(polskel)
        sizes[vid] = polskel.getLength()
        radii[vid] = lradii
    
    new_mtg.properties()['skeleton'] = skeletons
    new_mtg.properties()['skeletonradii'] = radii
    new_mtg.properties()['bbox'] = bbox
    new_mtg.properties()['size'] = sizes
    
    return new_mtg


from openalea.tree_matching.bipartitematching import *


def divide_range(vmin,vmax,subdiv, overlap = None):
    vrange = vmax - vmin
    vsubdivrange = vrange / subdiv
    if overlap is None:
        vsubdivrangeborder = vsubdivrange * 0.1
    else : vsubdivrangeborder = overlap
    vrangeset = [(vmin,vmin+vsubdivrange+vsubdivrangeborder)]
    vrangeset += [(vmin+i*vsubdivrange-vsubdivrangeborder, vmin+(i+1)*vsubdivrange+vsubdivrangeborder) for i in range(1,subdiv-1)]
    vrangeset += [(vmax-vsubdivrange-vsubdivrangeborder,vmax)]
    return vrangeset

def divide_bbox(bbx,subdiv = 6,overlap = None):
    xrangeset = divide_range(bbx.getXMin(),bbx.getXMax(),subdiv,overlap)
    yrangeset = divide_range(bbx.getYMin(),bbx.getYMax(),subdiv,overlap)
    zrangeset = divide_range(bbx.getZMin(),bbx.getZMax(),subdiv,overlap)
    bbxset = []
    for xmin,xmax in xrangeset:
        for ymin,ymax in yrangeset:
            for zmin,zmax in zrangeset:
                bbxset.append(BoundingBox((xmin,ymin,zmin),(xmax,ymax,zmax)))
    return bbxset


def create_edges_with_space_division(bbxprop1,bbxprop2,skeletonprop1,skeletonprop2,maxdistance):
    # Compute global bounding box
    globalbbox = compute_global_bbx(bbxprop2)
    
    #create space division to reach bbox
    bbxset = divide_bbox(globalbbox,subdiv=6,overlap=maxdistance)
    bbxcontains = {}
    for key,bbxkey in enumerate(bbxset):
        bbxcontains[key] = []
    for vid2, bbx2 in bbxprop2.iteritems():
        for key,bbxkey in enumerate(bbxset):
            if bbxkey.intersect(bbx2):
                bbxcontains[key] += [vid2]
                
        
    edges = []    
    for vid1, bbx in bbxprop1.iteritems():
        vid2tocompare = set()
        for key,bbxkey in enumerate(bbxset):
            if bbxkey.intersect(bbx):
                vid2tocompare = vid2tocompare.union(bbxcontains[key])
        for vid2 in vid2tocompare:
            if bbx.distance(bbxprop2[vid2]) <= maxdistance:
                edges.append( (vid1,vid2, skeletonprop1[vid1].hausdorff_distance(skeletonprop2[vid2])) )
    return edges

def create_edges(bbxprop1,bbxprop2,skeletonprop1,skeletonprop2,maxdistance):
    edges = []    
    for vid1, bbx in bbxprop1.iteritems():
        for vid2, bbx2 in bbxprop2.iteritems():
            if bbx.distance(bbx2) <= maxdistance:
                edges.append( (vid1,vid2, skeletonprop1[vid1].hausdorff_distance(skeletonprop2[vid2])) )
    return edges

def check_minimum_edges_cost(edges,set1,set2,delcost1,delcost2):
    me1, me2 = {}, {}
    delcost1 = dict([(i,c) for i,c in zip(set1,delcost1)])
    delcost2 = dict([(i,c) for i,c in zip(set2,delcost2)])
    
    for i,j,c in edges:
        me1[i] = min(me1.get(i,1e6),c)
        me2[j] = min(me2.get(j,1e6),c)
    
    pb1 = []
    for i in set1:
        if not me1.has_key(i):
            pb1.append(i)
        elif me1[i] > delcost1[i]:
            pb1.append(i)
    print len(pb1),pb1
    pb2 = []
    for j in set2:
        if not me2.has_key(j):
            pb2.append(j)
        elif me2[j] > delcost2[j]:
            pb2.append(j)
    print len(pb2),pb2
    
def create_comparison(mtg1,mtg2,maxdistance = None, cacheprefix = None, cachepath = get_shared_data('comparison')):    
    set1 = mtg1.vertices(scale=2)
    set2 = mtg2.vertices(scale=2)
    # bbox properties
    bbxprop1 = mtg1.property('bbox')
    bbxprop2 = mtg2.property('bbox')
    # skeleton
    skeletonprop1 = mtg1.property('skeleton')
    skeletonprop2 = mtg2.property('skeleton')
    # size
    sizeprop1 = mtg1.property('size')
    sizeprop2 = mtg2.property('size')
    # max distance
    if maxdistance is None:
        # take second element to avoid appended trunk at the begining
        firstelem = mtg1.roots_iter(scale=2).next()
        secondelem = mtg1.children_iter(firstelem).next()
        skel1 = skeletonprop1[secondelem]
        skel1length = Polyline(skel1).getLength()
        nbsegments = (len(skel1)-1)
        maxdistance = 2 * skel1length / nbsegments
        print 'Use max distance :',maxdistance,(skel1length , nbsegments)
    
    # Compute edge
    from time import clock
    start = clock()
    
    print 'Create edges of the bipartite graph'
    # edges = create_edges(bbxprop1,bbxprop2,skeletonprop1,skeletonprop2,maxdistance)
    edges = create_edges_with_space_division(bbxprop1,bbxprop2,skeletonprop1,skeletonprop2,maxdistance)
    nbedges = len(edges)
    print 'Number of edge and Average degree:',nbedges,',',nbedges/len(bbxprop1)
    
    print "Time to compute edges is : ", clock() - start
    if cacheprefix:
        edge_cache_file = os.path.join(cachepath,cacheprefix+'_edges.pkl')
        writefile(edge_cache_file,edges)
        
    use_max_size = False
    if use_max_size:    
        # Compute global bounding box
        globalbbox = compute_global_bbx(bbxprop2)    
        maxsize = norm(globalbbox.getSize())
        print 'Use max size as non matching cost :', maxsize
        delset1cost =  [maxsize for i in set1]
        delset2cost =  [maxsize for i in set2]

    else:
        delset1cost =  [sizeprop1[i] for i in set1]
        delset2cost =  [sizeprop2[i] for i in set2]

    # return BipartiteMatching
    return BipartiteMatching(set1,set2,edges,delset1cost,delset2cost)


def find_maximal_paths(ni,mtg,nonmapping):
    path = [ni]
    np = mtg.parent(ni)
    while np in nonmapping:
        path = [np]+path
        np = mtg.parent(np)
    nc = [c for c in mtg.children(ni) if c in nonmapping]
    if len(nc) == 0:
        if len(path) == 1: return None
        else : return [path]
    else:
        paths = [path+[c] for c in nc]
    endedpaths = []
    while len(paths) > 0:
        npaths = []
        for p in paths:
            nc = [c for c in mtg.children(p[-1]) if c in nonmapping]
            if len(nc) == 0:
                endedpaths.append(p)
            else:
                npaths += [p+[c] for c in nc]
        paths = npaths
    return endedpaths

def find_all_paths(node,maxpaths):
    if maxpaths is None: return None
    ind = maxpaths[0].index(node)
    result = []
    for maxpath in maxpaths:
        length = len(maxpath)
        for j in xrange(ind,length):
            for i in xrange(0,ind+1):
                candidate = maxpath[i:j+1]
                if not candidate in result:
                    result.append(candidate)
    return result

def uniondist(n1set, mtg1,  n2, mtg2):
    skel1 = mtg1.property('skeleton')
    ptslist = [p for p in skel1[n1set[0]]]
    for n1 in  n1set[1:]:
        ptslist += [p for p in skel1[n1]][1:]
    ptslist2 = mtg2.property('skeleton')[n2]
    return ptslist2.hausdorff_distance(Point3Array(ptslist))

def extend_comparison(ref_mtg,tested_mtg, mapping, nonmapping1, nonmapping2, verbose = False):
    """ A post treatment to find group of nodes that match node of the other graph """
    newmapping = []
    nonmapping1 = set(nonmapping1)
    nonmapping2 = set(nonmapping2)
    for ni,nj in mapping:
        actualdistvalue = ref_mtg.property('skeleton')[ni].hausdorff_distance(tested_mtg.property('skeleton')[nj])
        newmindist1 = 2*actualdistvalue
        newmindist2 = 2*actualdistvalue
        paths1 = find_all_paths(ni,find_maximal_paths(ni,ref_mtg,nonmapping1))
        if not paths1 is None:
            a = [ (i,uniondist(p,ref_mtg, nj, tested_mtg)) for i,p in enumerate(paths1)]
            a.sort(cmp=lambda x,y : cmp(x[1],y[1]))
            newmindistpath1 = paths1[a[0][0]]
            newmindist1 = a[0][1]
        paths2 = find_all_paths(nj,find_maximal_paths(nj,tested_mtg,nonmapping2))
        if not paths2 is None:
            a = [ (i,uniondist(p, tested_mtg, ni, ref_mtg)) for i,p in enumerate(paths2)]
            a.sort(cmp=lambda x,y : cmp(x[1],y[1]))
            newmindistpath2 = paths2[a[0][0]]
            newmindist2 = a[0][1]
        choices = [newmindist1,newmindist2,actualdistvalue]
        choice = choices.index(min(choices))
        if choice == 2:
            newmapping += [(ni,nj)]
        elif choice == 0:
            if verbose: print 'Extend mapping from',(ni,nj),' to ',(newmindistpath1,nj)
            newmapping += [(newmindistpath1,nj)]
            for v in newmindistpath1:
                if v != ni:
                    nonmapping1.discard(v)
        elif choice == 1:
            if verbose: print 'Extend mapping from',(ni,nj),' to ',(ni,newmindistpath2)
            newmapping += [(ni,newmindistpath2)]
            for v in newmindistpath2:
                if v != nj:
                    nonmapping2.discard(v)
    return newmapping, nonmapping1, nonmapping2

def flatten(l):
    res = []
    for v in l :
        if isinstance(v,list): res += flatten(v)
        else : res.append(v)
    return res

def topological_comparison(ref_mtg,tested_mtg, mapping, nonmapping1, nonmapping2, verbose = False):
    edge_mapping = []
    edge_nonmapping1 = []
    edge_nonmapping2 = []
    mapping1, mapping2 = {}, {}
    
    for nis,njs in mapping:
        if isinstance(nis,int):
            mapping1[nis] = njs
        else:
            for ni in nis: mapping1[ni] = njs
        if isinstance(njs,int):
            mapping2[njs] = nis
        else:
            for nj in njs: mapping2[nj] = nis
        
    for nis,njs in mapping:
        if isinstance(nis,int):
            ni_parent =  ref_mtg.parent(nis) 
        else:
            ni_parent = [ref_mtg.parent(i) for i in nis ]
            ni_parent = [i for i in ni_parent if not i in nis]    
            assert len(ni_parent) == 1
            ni_parent= ni_parent[0]
        if isinstance(njs,int):
            nj_parent = tested_mtg.parent(njs)
        else:
            nj_parent = [tested_mtg.parent(i) for i in njs ]
            nj_parent = [i for i in nj_parent if not i in njs]    
            assert len(nj_parent) == 1
            nj_parent= nj_parent[0]
        if ni_parent is None or nj_parent is None:
            if ni_parent == nj_parent: edge_mapping.append((nis,njs))
            else :
                edge_nonmapping1.append(nis)
                edge_nonmapping2.append(njs)
            continue
        while ni_parent in nonmapping1 : ni_parent = ref_mtg.parent(ni_parent)
        while nj_parent in nonmapping2 : nj_parent = tested_mtg.parent(nj_parent) 
        mappedniparent = mapping1.get(ni_parent,[]) 
        mappednjparent = mapping2.get(nj_parent,[]) 
        def is_same(pid,mappedid):
            if isinstance(mappedid,int): return mappedid == pid
            else : return pid in mappedid
            
        if is_same(nj_parent,mappedniparent):
            assert is_same(ni_parent,mappednjparent)
            edge_mapping.append((nis,njs))
        else:
            edge_nonmapping1.append(nis)
            edge_nonmapping2.append(njs)
    return edge_mapping, edge_nonmapping1,edge_nonmapping2
    

        
def comparison_process(ref_mtg,tested_mtg,cacheprefix = None,cachepath = get_shared_data('comparison')):
    """ Return score of geometrical comparison (nb of common elements) and topological comparison (nb of common edges) """
    from time import clock
    print 'mtgs sizes :',ref_mtg.nb_vertices(scale=ref_mtg.max_scale()),tested_mtg.nb_vertices(scale=tested_mtg.max_scale())
    print 'make mtgs comparable'
    ref_mtg = comparable_mtg(ref_mtg)
    tested_mtg = comparable_mtg(tested_mtg)
    print 'mtgs sizes :',ref_mtg.nb_vertices(scale=2),tested_mtg.nb_vertices(scale=2)
    
    if cacheprefix:
        if not os.path.exists(cachepath): os.makedirs(cachepath)    
        raw_comparison_cache_file = os.path.join(cachepath,cacheprefix+'_raw_comparison.pkl')
    
    if cacheprefix and os.path.exists(raw_comparison_cache_file) : 
        print 'Read raw matching from',repr(raw_comparison_cache_file)
        raw_comparison = readfile(raw_comparison_cache_file)
    else:
        print 'Compute matching'
        start = clock()
        matcher = create_comparison(ref_mtg,tested_mtg,None,cacheprefix,cachepath)
        if cacheprefix : 
            raw_comparison_cache_file_intermediate = os.path.join(cachepath,cacheprefix+'_raw_comparison_intermediate.pkl')
        else : raw_comparison_cache_file_intermediate = None
        raw_comparison = matcher.match(raw_comparison_cache_file_intermediate)
        print ">>>>>>>> The time for computing matching is : ", clock() - start , " <<<<<<<<. "
        if cacheprefix : writefile(raw_comparison_cache_file,raw_comparison)
    distance, mapping, nonmapping1, nonmapping2 = raw_comparison
    print 'Nb of mapping :',len(mapping), len(nonmapping1), len(nonmapping2)
    print 'Distance :',distance
    nbtotvertex = ref_mtg.nb_vertices(scale=2)+tested_mtg.nb_vertices(scale=2)
    nbmapped = 2.*len(mapping)
    score = nbmapped/float(nbtotvertex)
    print 'Mapping score :',score
    ref_sizes = ref_mtg.property('size')
    tested_sizes = tested_mtg.property('size')
    totlength = sum(ref_sizes.itervalues())+sum(tested_sizes.itervalues())
    matchedlength = sum([ref_sizes[i] for i,j in mapping])+sum([tested_sizes[j] for i,j in mapping])
    cscore = matchedlength/totlength
    print 'Weighted mapping score :',cscore
    
    # Extended comparison
    with_extension = False
    if with_extension:
        ext_comparison = extend_comparison(ref_mtg,tested_mtg,mapping, nonmapping1, nonmapping2)
        newmapping, newnonmapping1, newnonmapping2 = ext_comparison
        assert len(newmapping) == len(mapping)
        nbnewmapped = nbtotvertex-len(newnonmapping1)-len(newnonmapping2)
        print 'Nb element added to the mapping :',nbnewmapped-nbmapped
        nscore = nbnewmapped/float(nbtotvertex)
        print 'New mapping score :',nscore
        nweunmatchedlength = sum([ref_sizes[i] for i in newnonmapping1])+sum([tested_sizes[i] for i in newnonmapping2])
        newcscore = (totlength-nweunmatchedlength)/totlength
        print 'New Corrected mapping score :',newcscore
        if cacheprefix:
            ext_comparison_cache_file = os.path.join(cachepath,cacheprefix+'_extended_comparison.pkl')
            writefile(ext_comparison_cache_file,ext_comparison)      
        # Second visualization
        # # scene = comparison_representation(ref_mtg,tested_mtg,newmapping)
        # scene = ext_comparison_representation(ref_mtg,tested_mtg,mapping,newmapping)
        # Viewer.display(scene)
    else:
        newmapping, newnonmapping1, newnonmapping2 = mapping, nonmapping1, nonmapping2
        newcscore = cscore
    
    # Topological comparison
    edgematching ,nonedgematching1,nonedgematching2 = topological_comparison(ref_mtg,tested_mtg,newmapping, newnonmapping1, newnonmapping2)
    print map(len,[edgematching,nonedgematching1,nonedgematching2])
    tscore = (2*len(edgematching))/float(sum(map(len,[edgematching,edgematching,nonedgematching1,nonedgematching2])))
    print 'Topological comparison score',tscore
    return newcscore, tscore
    

