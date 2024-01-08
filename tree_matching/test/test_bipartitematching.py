from openalea.tree_matching.bipartitematching import BipartiteMatching
import time
from random import randint, uniform

def test(maxinrange = 3 ,maxoutrange = 2,deg = 1):
    print("test")
    set1 = [i for i in range(maxinrange)]
    set2 = [i+3 for i in range(maxoutrange)]
    edges = []
    for ni in set1:
        prevconnect = []
        for d in range(deg):
            newid = set2[randint(0,maxoutrange-1)]
            if not newid in prevconnect: # check for not connecting 2 times the same elements
                edges.append((ni,newid,uniform(0,10)))
                prevconnect.append(newid)
    print(edges)
    print((set1,set2))
    m = BipartiteMatching(set1,set2,edges,[20 for i in range(maxinrange)],[20 for i in range(maxoutrange)])
    print("Debut Matching") 
    start =  time.perf_counter()
    res = m.match()
    print(res)
    end =  time.perf_counter()
    print(("Elaspsed Time = ",end-start,"seconds"))
    print("Matching")
   # print res

def test2():
    maxinrange = 2 
    maxoutrange = 4
    set1 = [i+1 for i in range(maxinrange)]
    set2 = [i+4 for i in range(maxoutrange)]
    edges = [(1,4,4),(1,5,2),(2,5,1),(2,6,2),(1,7,10)]
    m = BipartiteMatching(set1,set2,edges,[5 for i in range(maxinrange)],[5 for i in range(maxoutrange)])
    start =  time.perf_counter()
    res = m.match()
    end =  time.perf_counter()
    print(res)
    assert res[1] == [(1, 5), (2, 6)]
    print(("Elaspsed Time = ",end-start,"seconds"))


def test3():
    maxinrange = 2 
    maxoutrange = 4
    set1 = [i+1 for i in range(maxinrange)]
    set2 = [i+4 for i in range(maxoutrange)]
    edges = [(1,4,3),(1,5,4),(2,4,1),(2,6,3),(2,7,10)]
    m = BipartiteMatching(set1,set2,edges,[5 for i in range(maxinrange)],[5 for i in range(maxoutrange)])
    start =  time.perf_counter()
    res = m.match()
    end =  time.perf_counter()
    print(res)
    assert res[1] == [(1, 5), (2, 4)]
    print(("Elaspsed Time = ",end-start,"seconds"))


def test_with_capacities():
    maxinrange = 2 
    maxoutrange = 4
    set1 = [i+1 for i in range(maxinrange)]
    set2 = [i+4 for i in range(maxoutrange)]
    edges = [(1,4,4),(1,5,2),(2,5,1),(2,6,2),(1,7,10)]
    capacities = [2,1]
    m = BipartiteMatching(set1,set2,edges,[10 for i in range(maxinrange)],[10 for i in range(maxoutrange)], capacities)
    start =  time.perf_counter()
    res = m.match()
    end =  time.perf_counter()
    print(res)
    assert res[1] == [(1, 4), (1, 5), (2, 6)]
    print(("Elaspsed Time = ",end-start,"seconds"))




if __name__ == '__main__':
    test()
    #test2()
    #test3()
    #test_with_capacities()