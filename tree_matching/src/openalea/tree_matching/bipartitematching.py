from .__tree_matching__ import GeneralMatchPath

def test_nb_theta_as_min(edges):
    a, b = {}, {}
    for ij,c in edges.items():
        i,j = ij
        if i != -1 and i in a:
            if a[i][0] > c: 
                a[i] = (c,j)
        else: a[i] = (c,j)
        if j != -1 and j in b:
            if b[j][0] > c: 
                b[j] = (c,i)
        else: b[j] = (c,i)
    print('Check', len([1 for c,j in a.values() if j == -1]),len([1 for c,i in b.values() if i == -1]))

def remove_elements(todelete,elements,nonmatchingcost):
    elements = list(elements)
    nonmatchingcost = list(nonmatchingcost)
    todelete.sort()
    todelete.reverse()
    for i in todelete:
        del elements[i]
        del nonmatchingcost[i]
    return elements, nonmatchingcost
    
def filter_orphans(set1,set2, edges, nonmatchingset1cost, nonmatchingset2cost):
    idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
    idmapset2 = dict([(nj,i) for i,nj in enumerate(set2)])
    
    degset1 = dict([(ni,0) for ni in set1])
    degset2 = dict([(nj,0) for nj in set2])
    
    for ni,nj,cost in edges:
        degset1[ni] += 1
        degset2[nj] += 1
            
    orphans1 = []
    todelete = []
    orphanscost = 0
    for ni,deg in degset1.items():
        if deg == 0 : 
            orphans1.append(ni)
            todelete.append(idmapset1[ni])
            orphanscost += nonmatchingset1cost[idmapset1[ni]]
    set1,nonmatchingset1cost = remove_elements(todelete, set1,nonmatchingset1cost)
            
    orphans2 = []
    todelete = []
    for nj,deg in degset2.items():
        if deg == 0 : 
            orphans2.append(nj)
            todelete.append(idmapset2[nj])
            orphanscost += nonmatchingset2cost[idmapset2[nj]]
            
    set2,nonmatchingset2cost = remove_elements(todelete, set2,nonmatchingset2cost)
    
    return set1,set2, nonmatchingset1cost, nonmatchingset2cost, orphans1, orphans2, orphanscost

def filter_insdel_sup_edges(set1,set2, edges, nonmatchingset1cost, nonmatchingset2cost):
    """ Filter edges with a cost sup than cost of deletion and insertion of its 2 nodes """
    delcost1 = dict(list(zip(set1,nonmatchingset1cost)))
    delcost2 = dict(list(zip(set2,nonmatchingset2cost)))
    
    newedges = []
    for ni,nj,cost in edges:
        insdelcost = delcost1[ni] + delcost2[nj]
        if insdelcost > cost:
            newedges += [(ni,nj,cost)]
    return newedges
    
def filter_perfect_matches(set1,set2, edges, nonmatchingset1cost, nonmatchingset2cost):
    """ Filter nodes that matches uniquely with a cost of 0 to a corresponding nodes """
    perfectly_matched1 = {}
    perfectly_matched2 = {}
    for ni,nj,cost in edges:
        if cost == 0:
            perfectly_matched1[ni] = perfectly_matched1.get(ni,[])+[nj]
            perfectly_matched2[nj] = perfectly_matched2.get(nj,[])+[ni]
    
    perfect_matches = []
    for ni,njs in perfectly_matched1.items():
        if len(njs) == 1 and perfectly_matched2[njs[0]] == [ni]: # A unique one-to-one relation
            perfect_matches.append((ni,njs[0]))
    
    if len(perfect_matches) == 0:
        return set1,set2, edges, nonmatchingset1cost, nonmatchingset2cost, perfect_matches
    
    toremove1 = set([i for i,j in perfect_matches])
    toremove2 = set([j for i,j in perfect_matches])

    newedges = []
    for ni,nj,cost in edges:
        if not ni in toremove1 and not nj in toremove2:
            newedges.append((ni,nj,cost))

    idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
    idmapset2 = dict([(nj,i) for i,nj in enumerate(set2)])
    
    toremove1id = [idmapset1[i] for i in toremove1]
    toremove2id = [idmapset2[i] for i in toremove2]
    
    set1,nonmatchingset1cost = remove_elements(toremove1id, set1,nonmatchingset1cost)
    set2,nonmatchingset2cost = remove_elements(toremove2id, set2,nonmatchingset2cost)
             
    return set1,set2, newedges, nonmatchingset1cost, nonmatchingset2cost, perfect_matches

    
class BipartiteMatching:
    def __init__(self, set1, set2, possiblematching, nonmatchingset1cost, nonmatchingset2cost, set1capacity = None, filter = True):
        """ Built with the ids of the two sets and the possible matching between them (list of triplet (id1,id2,cost).
            nonmatchingset1cost is an ordered list containing the costs to not match an element of set 1. Same with set2 and nonmatchingset2cost.  """
                    
        assert len(nonmatchingset1cost) == len(set1) and 'One non matching cost for one element of set1'
        assert len(nonmatchingset2cost) == len(set2) and 'One non matching cost for one element of set2'
        
        self.orphans1, self.orphans2, self.perfect_matches = [], [], []
        self.orphanscost = 0
        if filter:
            ipossiblematching = possiblematching
            possiblematching = filter_insdel_sup_edges(set1, set2, ipossiblematching, nonmatchingset1cost, nonmatchingset2cost)
            # print len(ipossiblematching)-len(possiblematching)
            iipossiblematching = possiblematching
            set1,set2, possiblematching, nonmatchingset1cost, nonmatchingset2cost, self.perfect_matches = filter_perfect_matches(set1, set2, iipossiblematching, nonmatchingset1cost, nonmatchingset2cost)
            # print len(self.perfect_matches), len(iipossiblematching)-len(possiblematching)
            set1,set2, nonmatchingset1cost, nonmatchingset2cost, self.orphans1, self.orphans2, self.orphanscost = filter_orphans(set1, set2, possiblematching, nonmatchingset1cost, nonmatchingset2cost)
            # print len(self.orphans1), len(self.orphans2)
            # print 'Orphans :',self.orphans1 , self.orphans2
        
        self.nb_ni = len(set1)
        self.nb_nj = len(set2)
        
        assert len(nonmatchingset1cost) == self.nb_ni and 'Require one non matching cost for each element of set1'
        assert len(nonmatchingset2cost) == self.nb_nj and 'Require one non matching cost for each element of set2'
        
        # map from set1 id to normalized id
        idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
        idmapset2 = dict([(nj,self.nb_ni+i) for i,nj in enumerate(set2)])
        
        # reverse map to use to translate back results
        self.idmapset = dict([(i,ni) for ni,i in idmapset1.items()]+[(i,nj) for nj,i in idmapset2.items()])
        self.idmapset[-1] = None
        
        edges_cost = {}
        inconnect = [[] for it in range(self.nb_ni)]
        outconnect = [[] for it in range(self.nb_nj)]
        for ni,nj,cost in possiblematching:
            # normalize index
            nni, nnj = idmapset1[ni],idmapset2[nj]
            # fill in and out connections
            inconnect[nni].append( nnj)
            outconnect[nnj-self.nb_ni].append( nni)
            # create normalized edge cost map
            edges_cost[(nni,nnj)] = cost
            
        # insert cost of non matching        
        for ni, cost in zip(set1,nonmatchingset1cost):
            edges_cost[(idmapset1[ni],-1)] = cost
        for nj, cost in zip(set2,nonmatchingset2cost):
            edges_cost[(-1,idmapset2[nj])] = cost
        
        edges_cost[(-1,-1)] = 0
        self.edges_cost = edges_cost
        #test_nb_theta_as_min(edges_cost)
        
        # Initialize the Cpp structure
        class GeneralMatchPathWrapper(GeneralMatchPath):
            def __init__(self,ni,nj,inconnect,outconnect,edges_cost):
                if set1capacity:
                    GeneralMatchPath.__init__(self,list(range(ni)),list(range(nj)),inconnect,outconnect, set1capacity)
                else:
                    GeneralMatchPath.__init__(self,list(range(ni)),list(range(nj)),inconnect,outconnect)
                self.edges_cost = edges_cost
               
            def edgeCost(self,a,b):
                return self.edges_cost[(a,b)]

        self.__cppmatchpath = GeneralMatchPathWrapper(self.nb_ni,self.nb_nj,inconnect,outconnect,edges_cost)
        
    def match(self,cachefile = None):
        distance, matching = self.__cppmatchpath.bipartiteMatching()
        # print 'Nb matching found :',len(matching) # Attention ceci compte les paires "non matchees",
                                                  # certaines sont comptees deux fois
        if cachefile:
            import pickle as pickle
            pickle.dump(matching, file(cachefile,'wb'), pickle.HIGHEST_PROTOCOL)
        newmatching = set()
        nonmatching1 = set()
        nonmatching2 = set()
        for ni,nj in matching:
            #[Pascal] J'ai corrige pour ne voir apparaitre qu'une seule fois un noeud non mappe
            if ni == -1 or nj == -1:
               if ni == -1: ni = nj   # put value diff from -1 in ni
               if ni == -1: continue  # if both are -1, do not consider the tuple
               if ni >= self.nb_ni: nonmatching2.add(self.idmapset[ni]) # if ni > nb_ni, this means that it is part of set2
               else : nonmatching1.add(self.idmapset[ni])  # else it is part of set1
            else:
                if ni >= self.nb_ni: ni,nj = nj,ni
                newelem = (self.idmapset[ni],self.idmapset[nj])                
                newmatching.add(newelem)
        # print             

        return distance+self.orphanscost, list(newmatching)+self.perfect_matches , self.orphans1+list(nonmatching1), self.orphans2+list(nonmatching2)