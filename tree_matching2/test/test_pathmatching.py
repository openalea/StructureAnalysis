from openalea.tree_matching import *


class MyMatching(GeneralMatchPath):
    def __init__(self,maxinrange,maxoutrange,deg):
        inconnect = [range(maxinrange+min(maxoutrange,max(0,i-deg/2)),maxinrange+min(maxoutrange,1+i+deg/2)) for i in xrange(maxinrange)]
        outconnect =  [range(min(maxinrange,max(0,i-deg/2)), min(maxinrange,1+i+deg/2)) for i in xrange(maxoutrange)]
        GeneralMatchPath.__init__(self,range(maxinrange),range(maxinrange,maxinrange+maxoutrange),inconnect,outconnect)
        print inconnect
        print outconnect
    def edgeCost(self,a,b):
            print a,b
            return abs(a-b)


def test(maxinrange = 15 ,maxoutrange = 10,deg = 3):
    m = MyMatching(maxinrange,maxoutrange,deg)
    res = m.bipartiteMatching()
    print res

if __name__ == '__main__':
    test()