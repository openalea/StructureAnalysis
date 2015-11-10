from openalea.tree_matching import *
import time

class MyMatching(GeneralMatchPath):
    def __init__(self,maxinrange,maxoutrange,deg):
        inconnect = [range(maxinrange+min(maxoutrange,max(0,i-deg/2)),maxinrange+min(maxoutrange,1+i+deg/2)) for i in xrange(maxinrange)]
        outconnect =  [range(min(maxinrange,max(0,i-deg/2)), min(maxinrange,1+i+deg/2)) for i in xrange(maxoutrange)]
        GeneralMatchPath.__init__(self,range(maxinrange),range(maxinrange,maxinrange+maxoutrange),inconnect,outconnect)
        #print "inconnect = " , inconnect
        #print "outconnect = " , outconnect
    def edgeCost(self,a,b):
            #print "EdgeCost = ",a,b
            if a == -1 or b == -1:
                return 0
            else:
                return abs(a-(b-100))


def test(maxinrange = 100 ,maxoutrange = 100,deg = 3):
    print "test"
    m = MyMatching(maxinrange,maxoutrange,deg)
    print "Debut Matching" 
    start = time.clock()
    res = m.bipartiteMatching()
    print res
    end = time.clock()
    print "Elaspsed Time = ",end-start,"seconds"
    print "Matching"
   # print res

if __name__ == '__main__':
    test()
