from __tree_matching__ import GeneralMatchPath

def test_nb_theta_as_min(edges):
    a, b = {}, {}
    for ij,c in edges.iteritems():
        i,j = ij
        if i != -1 and a.has_key(i):
            if a[i][0] > c: 
                a[i] = (c,j)
        else: a[i] = (c,j)
        if j != -1 and b.has_key(j):
            if b[j][0] > c: 
                b[j] = (c,i)
        else: b[j] = (c,i)
    print 'Check', len([1 for c,j in a.itervalues() if j == -1]),len([1 for c,i in b.itervalues() if i == -1])

def filter_elements(set1,set2, edges, nonmatchingset1cost, nonmatchingset2cost):
    idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
    idmapset2 = dict([(nj,i) for i,nj in enumerate(set2)])
    degset1 = dict([(ni,0) for i,ni in enumerate(set1)])
    degset2 = dict([(nj,0) for i,nj in enumerate(set2)])
    newedges = []
    
    for ni,nj,cost in edges:
        insdelcost = nonmatchingset1cost[idmapset1[ni]] + nonmatchingset2cost[idmapset2[nj]]
        if insdelcost > cost:
            newedges += [(ni,nj,cost)]
            degset1[ni] += 1
            degset2[nj] += 1
            
    nonmatching1 = []
    todelete = []
    for ni,deg in degset1.iteritems():
        if deg == 0: 
            nonmatching1.append(ni)
            todelete.append(idmapset1[ni])
    
    todelete.sort()
    todelete.reverse()
    for i in todelete:
        del set1[i]
        del nonmatchingset1cost[i]
            
    nonmatching2 = []
    todelete = []
    for nj,deg in degset2.iteritems():
        if deg == 0: 
            nonmatching2.append(nj)
            todelete.append(idmapset2[nj])

    todelete.sort()
    todelete.reverse()
    for i in todelete:
        del set2[i]
        del nonmatchingset2cost[i]
             
    return set1,set2, newedges, nonmatchingset1cost, nonmatchingset2cost, nonmatching1, nonmatching2

class BipartiteMatching:
    def __init__(self, set1, set2, ipossiblematching, nonmatchingset1cost, nonmatchingset2cost):
        """ Built with the ids of the two sets and the possible matching between them (list of triplet (id1,id2,cost).
            nonmatchingset1cost is an ordered list containing the costs to not match an element of set 1. Same with set2 and nonmatchingset2cost.  """
        assert len(nonmatchingset1cost) == len(set1) and 'One non matching cost for one element of set1'
        assert len(nonmatchingset2cost) == len(set2) and 'One non matching cost for one element of set2'
        
 #       set1,set2, possiblematching, nonmatchingset1cost, nonmatchingset2cost, self.orphans1, self.orphans2 = filter_elements(set1, set2, ipossiblematching, nonmatchingset1cost, nonmatchingset2cost)
 #       print len(self.orphans1), len(self.orphans2), len(ipossiblematching)-len(possiblematching)
        
        self.nb_ni = len(set1)
        self.nb_nj = len(set2)
        
        # map from set1 id to normalized id
        idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
        idmapset2 = dict([(nj,self.nb_ni+i) for i,nj in enumerate(set2)])
        
        # reverse map to use to translate back results
        self.idmapset = dict([(i,ni) for ni,i in idmapset1.iteritems()]+[(i,nj) for nj,i in idmapset2.iteritems()])
        self.idmapset[-1] = None
        
        edges_cost = {}
        inconnect = [[] for it in xrange(self.nb_ni)]
        outconnect = [[] for it in xrange(self.nb_nj)]
  #      for ni,nj,cost in possiblematching:
        for ni,nj,cost in ipossiblematching:
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
        
        #test_nb_theta_as_min(edges_cost)
        
        # Initialize the Cpp structure
        class GeneralMatchPathWrapper(GeneralMatchPath):
            def __init__(self,ni,nj,inconnect,outconnect,edges_cost):
                GeneralMatchPath.__init__(self,range(ni),range(nj),inconnect,outconnect)
                self.edges_cost = edges_cost
               
            def edgeCost(self,a,b):
                #print a,b,self.edges_cost[(a,b)]
                return self.edges_cost[(a,b)]

        self.__cppmatchpath = GeneralMatchPathWrapper(self.nb_ni,self.nb_nj,inconnect,outconnect,edges_cost)
        
    def match(self,cachefile = None):
        distance, matching = self.__cppmatchpath.bipartiteMatching()
        print 'Nb matching found :',len(matching) # Attention ceci compte les paires "non matchees",
                                                  # certaines sont comptees deux fois
        if cachefile:
            import cPickle as pickle
            pickle.dump(matching, file(cachefile,'wb'), pickle.HIGHEST_PROTOCOL)
        newmatching = set()
        nonmatching1 = [ ]
        nonmatching2 = [ ]
        for ni,nj in matching:
            #print ni,nj,'  ',
            #[Pascal] J'ai corrige pour ne voir apparaitre qu'une seule fois un noeud non mappe
            if ni == -1 or nj == -1:
                if ni != -1 and ni < self.nb_ni:
                    nonmatching1.append(self.idmapset[ni])
                else:
                    if nj != -1 and nj >= self.nb_ni:
                        nonmatching2.append(self.idmapset[nj])
#                if ni == -1: ni,nj = nj,ni
#                if ni >= self.nb_ni: nonmatching2.append(self.idmapset[ni])
#                else : nonmatching1.append(self.idmapset[ni])
            else:
                if ni >= self.nb_ni: ni,nj = nj,ni
                newelem = (self.idmapset[ni],self.idmapset[nj])                
                newmatching.add(newelem)
                     

        return distance, list(newmatching) ,nonmatching1, nonmatching2 #, self.orphans1+nonmatching1, self.orphans2+nonmatching2
