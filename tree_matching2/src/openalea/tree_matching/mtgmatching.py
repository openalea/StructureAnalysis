from openalea.tree_matching import Matching, ExtMatching, NodeCost
from mtgimport import mtg2treegraph


class NodeCostWrapper (NodeCost):
    def __init__(self, mtgnodecost, idmap1, idmap2) :
        NodeCost.__init__(self)
        self.idmap1 = idmap1
        self.idmap2 = idmap2
        self.mtgnodecost = mtgnodecost
        
    def getDeletionCost(self,a) :
        return self.mtgnodecost.getDeletionCost(self.idmap1[a.id])
        
    def getInsertionCost(self,b) :
        return self.mtgnodecost.getInsertionCost(self.idmap2[b.id])
      
    def getChangingCost(self,a,b) :
        return self.mtgnodecost.getChangingCost(self.idmap1[a.id],self.idmap2[b.id])
     
    def getMergingCost(self,a,b) :
        return self.mtgnodecost.getMergingCost([self.idmap1[i.id] for i in a],self.idmap2[b.id])

    def getSplittingCost(self,a,b) :
        return self.mtgnodecost.getSplittingCost(self.idmap1[a.id],[self.idmap2[i.id] for i in b])


class MtgMatchingWrapper :
    def __init__(self,mtg1,mtg2,scale1,scale2,cost,tabletype = 0,root1=None,root2=None):
        # mtg value
        self.mtg1 = mtg1
        self.mtg2 = mtg2
        # creation of coresponding tree graph from mtg
        self.tree1,self.idmap1 = mtg2treegraph(self.mtg1,scale1,root1)
        self.tree2,self.idmap2 = mtg2treegraph(self.mtg2,scale2,root2)
        # id map from treegraph to mtg 
        self.idmap1 = dict([(j,i) for i,j in self.idmap1.iteritems()])
        self.idmap2 = dict([(j,i) for i,j in self.idmap2.iteritems()])
        # making cost function.
        self.cost = cost
        self.costwrapper = NodeCostWrapper(self.cost,self.idmap1,self.idmap2)
    
    def __getroot(self,a,b):
        # get root id
        if a is None : a = 0
        else : a = self.idmap1[a]
        if b is None : b = 0
        else : b = self.idmap2[b]
        return a,b
        
    def getList(self,a = None, b = None):
        a,b = self.__getroot(a,b) 
        res = self.__native_get_list(a,b)
        return [(self.idmap1[i],self.idmap2[j],k) for i,j,k in res]
    
    def __native_get_list(self,a,b):
        # to be redefined to point toward actual getList of cpp implementation
        return Matching.getList(self,a,b)
    
    def getDBT(self,a = None, b = None):
        a,b = self.__getroot(a,b)
        return self.__native_get_dbt(a,b)

    def __native_get_dbt(self,a,b):
        # to be redefined to point toward actual getDBT of cpp implementation
        return Matching.getDBT(self,a,b)

class MtgMatching (MtgMatchingWrapper,Matching):

    def __init__(self,mtg1,mtg2,scale1,scale2,cost,tabletype = 0,root1=None,root2=None):
        # mtg wrapper
        MtgMatchingWrapper.__init__(self,mtg1,mtg2,scale1,scale2,cost,tabletype,root1,root2)
        # creation of Matching
        Matching.__init__(self,self.tree1,self.tree2,self.costwrapper,tabletype)

    def __native_get_list(self,a,b):
        # to be redefined to point toward actual getList of cpp implementation
        return Matching.getList(self,a,b)

    def __native_get_dbt(self,a,b):
        # to be redefined to point toward actual getDBT of cpp implementation
        return Matching.getDBT(self,a,b)

class MtgExtMatching (MtgMatchingWrapper,ExtMatching):

    def __init__(self,mtg1,mtg2,scale1,scale2,cost,tabletype = 0,root1=None,root2=None):
        # mtg wrapper
        MtgMatchingWrapper.__init__(self,mtg1,mtg2,scale1,scale2,cost,tabletype,root1,root2)
        # creation of Matching
        ExtMatching.__init__(self,self.tree1,self.tree2,self.costwrapper,tabletype)

    def __native_get_list(self,a,b):
        # to be redefined to point toward actual getList of cpp implementation
        return ExtMatching.getList(self,a,b)

    def __native_get_dbt(self,a,b):
        # to be redefined to point toward actual getDBT of cpp implementation
        return ExtMatching.getDBT(self,a,b)

