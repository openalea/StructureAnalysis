from openalea.tree_matching import *
import time

class MyMatching(GeneralMatchPath):
    def __init__(self,maxinrange,maxoutrange,deg):
        inconnect = [list(range(int(maxinrange+min(maxoutrange,max(0,i-deg/2))), 
                                int(maxinrange+min(maxoutrange,1+i+deg/2)) ))  for i in range(maxinrange)]
        outconnect =  [list(range(int(min(maxinrange,max(0,i-deg/2))), 
                                  int(min(maxinrange,1+i+deg/2)) )) for i in range(maxoutrange)]
        GeneralMatchPath.__init__(self,list(range(maxinrange)),list(range(maxinrange,maxinrange+maxoutrange)),inconnect,outconnect)
        print("inconnect = " , inconnect)
        print("outconnect = " , outconnect)
    def edgeCost(self,a,b):
            #print "EdgeCost = ",a,b
            if a == -1 or b == -1:
                return 1.
            else:
                return abs(a-b)


def test(maxinrange = 3 ,maxoutrange = 3,deg = 2):
    print("test")
    m = MyMatching(maxinrange,maxoutrange,deg)
    print("Debut Matching") 
    start = time.perf_counter()
    res = m.bipartiteMatching()
    end = time.perf_counter()
    print("Elaspsed Time = ",end-start,"seconds")
    print("Matching")
    print(res)

if __name__ == '__main__':
    test()
