from openalea.tree_matching.bipartitematching import *
import time
from random import randint, uniform

def test(maxinrange = 10 ,maxoutrange = 10,deg = 2):
    print "test"
    set1 = [randint(0,1e9) for i in xrange(maxinrange)]
    set2 = [randint(0,1e9) for i in xrange(maxoutrange)]
    edges = []
    for ni in set1:
        prevconnect = []
        for d in xrange(deg):
            newid = set2[randint(0,maxoutrange-1)]
            if not newid in prevconnect: # check for not connecting 2 times the same elements
                edges.append((ni,newid,uniform(0,10)))
                prevconnect.append(newid)
    print edges
    m = BipartiteMatching(set1,set2,edges,[20 for i in xrange(maxinrange)],[20 for i in xrange(maxoutrange)])
    print "Debut Matching" 
    start = time.clock()
    res = m.match()
    print res
    end = time.clock()
    print "Elaspsed Time = ",end-start,"seconds"
    print "Matching"
   # print res

if __name__ == '__main__':
    test()