from __tree_matching__ import GeneralMatchPath


class BipartiteMatching:
    def __init__(self, set1, set2, possiblematching, nonmatchingset1cost, nonmatchingset2cost):
        """ Built with the ids of the two sets and the possible matching between them (list of triplet (id1,id2,cost).
            nonmatchingset1cost is an ordered list containing the costs to not match an element of set 1. Same with set2 and nonmatchingset2cost.  """
            
        self.nb_ni = len(set1)
        self.nb_nj = len(set2)
        
        assert len(nonmatchingset1cost) == self.nb_ni and 'One non matching cost for one element of set1'
        assert len(nonmatchingset2cost) == self.nb_nj and 'One non matching cost for one element of set2'
        
        # map from set1 id to normalized id
        idmapset1 = dict([(ni,i) for i,ni in enumerate(set1)])
        idmapset2 = dict([(nj,self.nb_ni+i) for i,nj in enumerate(set2)])
        
        # reverse map to use to translate back results
        self.idmapset = dict([(i,ni) for ni,i in idmapset1.iteritems()]+[(i,nj) for nj,i in idmapset2.iteritems()])
        self.idmapset[-1] = None
        
        edges_cost = {}
        inconnect = [[] for it in xrange(self.nb_ni)]
        outconnect = [[] for it in xrange(self.nb_nj)]
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
            
        # Initialize the Cpp structure
        class GeneralMatchPathWrapper(GeneralMatchPath):
            def __init__(self,ni,nj,inconnect,outconnect,edges_cost):
                GeneralMatchPath.__init__(self,range(ni),range(nj),inconnect,outconnect)
                self.edges_cost = edges_cost
               
            def edgeCost(self,a,b):
                return self.edges_cost[(a,b)]

        self.__cppmatchpath = GeneralMatchPathWrapper(self.nb_ni,self.nb_nj,inconnect,outconnect,edges_cost)
        
    def match(self,cachefile = None):
        distance, matching = self.__cppmatchpath.bipartiteMatching()
        print 'Nb matching found :',len(matching)
        if cachefile:
            import cPickle as pickle
            pickle.dump(matching, file(cachefile,'wb'), pickle.HIGHEST_PROTOCOL)
        newmatching = set()
        nonmatching1 = [ ]
        nonmatching2 = [ ]
        for ni,nj in matching:
            if ni == -1 or nj == -1:
                if ni == -1: ni,nj = nj,ni
                if ni >= self.nb_ni: nonmatching2.append(self.idmapset[ni])
                else : nonmatching1.append(self.idmapset[ni])
            else:
                if ni >= self.nb_ni: ni,nj = nj,ni
                newelem = (self.idmapset[ni],self.idmapset[nj])                
                newmatching.add(newelem)
        return distance, list(newmatching), nonmatching1, nonmatching2