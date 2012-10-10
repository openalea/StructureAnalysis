# -*- coding: cp1252 -*-

#from amlPy import *

from vplants.self_similarity.graph import *
from vplants.self_similarity.cst import *
from vplants.self_similarity.mtg2graph import *
from vplants.tree_reduction.graph import *
from vplants.tree_reduction.functions import *
from openalea.container.graph import *
import string



"""   dichotomic tree
def test1():
    edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(10,11),(11,12),(12,13),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,20)]
    edges_weight={0:2,1:2,2:2,3:2,4:2,5:2,6:2,7:2,8:2,9:2,10:2,11:2,12:2,13:2,14:2,15:2,16:2,17:2,18:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 20
    assert g.nb_edges() == 19
    return g"""

""" multiscale plant"""
def test1():
    edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(6,12),(7,8),(8,9),(9,10),(10,11),(12,13),(13,14),(14,15),(15,8)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 15
    assert g.nb_edges() == 15
    return g


""" intersected paths of figure 7
def test1():
    edges=[(1,2),(1,4),(2,3),(2,4),(3,4),(3,6),(4,6),(4,5),(5,6)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 6
    assert g.nb_edges() == 9
    return g"""

""" linear paths
def test1():
    edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 9
    assert g.nb_edges() == 8
    return g"""

""" first example
def test1():
    edges=[(1,2),(1,3),(1,4),(2,3),(2,5),(3,4),(3,5),(4,5)]
    edges_weight={0:1,1:2,2:1,3:1,4:1,5:1,6:1,7:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 5
    assert g.nb_edges() == 8
    return g"""

""" perfectly included paths of figure 7
def test1():
    edges=[(1,2),(1,8),(2,3),(2,8),(3,4),(3,7),(4,5),(4,8),(5,6),(5,8),(6,7),(6,9),(7,9),(8,9)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:2,13:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 9
    assert g.nb_edges() == 14
    return g"""

""" perfectly included paths2
def test1():
    edges=[(1,2),(1,8),(2,3),(2,8),(3,4),(3,7),(4,5),(4,8),(5,6),(5,8),(6,7),(6,9),(7,9),(8,9),(10,1),(10,11),(11,2)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:2,13:1,14:1,15:1,16:1}
    root = 10
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 11
    assert g.nb_edges() == 17
    return g"""

"""   rice panicle
def test1():
    edges=[(1,2),(2,3),(2,4),(3,5),(4,6),(4,7),(5,8),(5,9),(6,10),(7,11),(7,12),(8,13),(9,20),(9,25),(10,106),(11,53),(12,64),(12,74),(13,14),(14,15),(15,16),(16,17),(17,18),(18,19),(19,100),(20,21),(20,28),(21,22),(22,23),(23,24),(24,94),(24,101),(25,26),(26,27),(27,28),(28,104),(29,30),(30,31),(31,32),(32,33),(33,34),(34,35),(35,93),(35,89),(36,37),(36,45),(37,38),(38,39),(39,40),(40,41),(41,42),(42,43),(42,101),(43,44),(43,93),(44,93),(45,46),(46,47),(47,48),(48,49),(49,50),(50,51),(51,52),(51,90),(52,89),(52,90),(53,28),(53,54),(54,55),(54,104),(55,56),(56,57),(57,58),(58,59),(59,60),(59,101),(60,61),(60,93),(61,62),(61,93),(62,63),(62,88),(63,87),(63,101),(64,65),(65,66),(65,28),(66,67),(67,68),(68,69),(69,70),(70,71),(70,93),(71,72),(71,99),(72,73),(72,99),(73,87),(73,88),(74,75),(75,76),(76,77),(77,78),(78,79),(79,80),(80,81),(81,82),(81,100),(82,83),(82,98),(83,84),(83,90),(84,85),(84,90),(85,86),(85,88),(86,101),(86,100),(87,101),(87,103),(88,93),(88,102),(89,90),(89,91),(90,88),(90,102),(91,92),(91,102),(93,101),(93,104),(94,95),(94,93),(95,96),(95,99),(96,97),(96,98),(97,98),(97,100),(98,100),(98,102),(99,102),(99,87),(100,93),(100,103),(101,102),(101,104),(102,103),(102,104),(103,104),(104,105)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1,15:1,16:1,17:1,18:1,19:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,30:1,31:1,32:1,33:1,34:1,35:1,36:1,37:1,38:1,39:1,40:1,41:1,42:1,43:1,44:1,45:1,46:1,47:1,48:1,49:1,50:1,51:1,52:1,53:1,54:1,55:1,56:1,57:1,58:1,59:1,60:1,61:1,62:1,63:1,64:1,65:1,66:1,67:1,68:1,69:1,70:1,71:1,72:1,73:1,74:1,75:1,76:1,77:1,78:1,79:1,80:1,81:1,82:1,83:1,84:1,85:1,86:1,87:1,88:1,89:1,90:1,91:1,92:1,93:1,94:1,95:1,96:1,97:1,98:1,99:1,100:1,101:1,102:1,103:1,104:1,105:1,106:1,107:1,108:1,109:1,110:1,111:1,112:1,113:1,114:1,115:1,116:1,117:1,118:1,119:1,120:1,121:1,122:1,123:1,124:1,125:1,126:1,127:1,128:1,129:1,130:1,131:1,132:1,133:1,134:1,135:1,136:1,137:1,138:1,139:1,140:1,141:1,142:1,143:1,144:1,145:1,146:1,147:1,148:1,149:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 106
    assert g.nb_edges() == 150
    return g"""


"""   NEST of rice panicle
def test1():
    edges=[(1,2),(2,3),(2,7),(3,4),(3,5),(4,5),(4,8),(5,6),(5,9),(6,7),(7,8),(7,13),(8,9),(8,11),(9,10),(9,20),(10,11),(10,23),(11,12),(12,13),(13,14),(13,19),(14,15),(14,18),(15,16),(15,18),(16,17),(16,18),(17,18),(18,19),(19,20),(20,21),(20,23),(21,22),(21,24),(22,23),(22,24),(23,24),(24,25)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1,15:1,16:1,17:1,18:1,19:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,30:1,31:1,32:1,33:1,34:1,35:1,36:1,37:1,38:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 25
    assert g.nb_edges() == 39
    return g"""

"""    T0
def test1():
    edges=[(1,2) ,(1,4) ,(2,3) ,(2,4) ,(3,4) ,(4,5),(4,7) ,(5,6) ,(5,7) ,(6,7) ,(7,8) ,(7,10) ,(8,9) ,(8,10) ,(9,10) ,(10,11)]
    edges_weight={0:1,1:1,2:1,3:1,4:2,5:1,6:1,7:1,8:1,9:2,10:1,11:1,12:1,13:1,14:2,15:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 11
    assert g.nb_edges() == 16
    return g"""

"""    T1
def test1():
    edges=[(1,2),(1,4),(2,3),(2,4),(3,4),(4,5),(4,7),(5,6),(5,7),(6,7),(7,8),(7,10),(8,9),(8,10),(9,10),(10,11),(10,12),(11,12),(12,13),(12,14),(13,14),(14,15),(14,16),(15,16),(16,17)]
    edges_weight={0:1,1:1,2:1,3:1,4:2,5:1,6:1,7:1,8:1,9:2,10:1,11:1,12:1,13:1,14:2,15:1,16:1,17:2,18:1,19:1,20:2,21:1,22:1,23:2,24:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 17
    assert g.nb_edges() == 25
    return g"""

"""    T2
def test1():
    edges=[(1,2),(1,10),(2,3),(2,10),(3,4),(3,10),(4,5),(4,7),(5,6),(5,7),(6,7),(7,8),(7,10),(8,9),(8,10),(9,10),(10,11),(10,12),(11,12),(12,13),(12,14),(13,14),(14,15),(14,16),(15,16),(16,17)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:2,11:1,12:1,13:1,14:1,15:2,16:1,17:1,18:2,19:1,20:1,21:2,22:1,23:1,24:2,25:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 17
    assert g.nb_edges() == 26
    return g"""

"""    T3
def test1():
    edges=[(1,2),(1,13),(2,3),(3,4),(3,13),(4,5),(5,6),(6,13),(6,7),(7,8),(7,13),(8,9),(9,10),(9,13),(10,11),(11,12),(11,13),(12,13)]
    edges_weight={0:1,1:1,2:3,3:1,4:1,5:3,6:1,7:1,8:3,9:1,10:1,11:3,12:1,13:1,14:3,15:1,16:1,17:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 13
    assert g.nb_edges() == 18
    return g"""

"""    palmier
def test1():
    edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(6,10),(7,8),(7,11),(8,9),(8,12),(9,13),(10,11),(11,12),(12,13)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:2,7:1,8:2,9:1,10:2,11:3,12:1,13:1,14:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 13
    assert g.nb_edges() == 15
    return g"""

"""    Borchet and Honda model
def test1():
    edges=[(1,2),(2,3),(2,14),(3,4),(3,14),(4,5),(4,14),(5,6),(5,14),(6,7),(6,14),(7,8),(7,14),(8,9),(8,14),(9,10),(9,14),(10,11),(10,14),(11,12),(11,14),(12,13),(12,15),(13,16),(14,15),(14,17),(15,16),(15,17),(16,17)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1,15:1,16:1,17:1,18:1,19:1,20:1,21:1,22:1,23:2,24:1,25:1,26:1,27:1,28:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 17
    assert g.nb_edges() == 29
    return g"""

""" arbre vide
def test1():
    edges=[]
    edges_weight={}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 0
    assert g.nb_edges() == 0
    return g"""


def taille_arbre(g,vid):
    if descendents(vid, g)==[]:
        return 1
    else:
        taillev=1
        for v2 in g.out_neighbors(vid):
            eid=find_eid(vid, v2, g)
            poid=g._edge_property['weight'][eid]
            taillev=taillev+(poid*taille_arbre(g,v2))
        return taillev
            

def mapping_cgraph_test1(graph):
    """
    transformer la structure du DAG de la classe c_graph vers la structure test1"""
    
    edges=[]
    edges_weight={}
    k=0
    for i in range (0, nb_nodes(graph)):
        n = node(graph, i)
        """print "Edges From Node ", n.label
        print "Signature",n.signature"""
        for j in range(nb_edges(n)):
            e = edge(n, j)
            if e.start == n:
                edges=[(nb_nodes(graph)-n.signature,nb_nodes(graph)-e.end.signature)]+edges
                edges_weight[k]=int(e.label)
                k=k+1
                """print  n.label, " ---", e.label,"--->", e.end.label"""

    """
    for i in range(len(graph.edges)):
        print [(graph.nb_signature-graph.edges[i].start.signature,graph.nb_signature-graph.edges[i].end.signature)]
        edges=[(graph.nb_signature-graph.edges[i].start.signature,graph.nb_signature-graph.edges[i].end.signature)]+edges
        edges_weight[i]=int(graph.edges[i].label)"""
    root=1        
    g = from_edges(root, edges, edges_weight)
    print nb_nodes(graph), k
    if (nb_nodes(graph)>1):
        assert g.nb_vertices() == nb_nodes(graph)
        assert g.nb_edges() == k
        return g
    else:
        return g
        

def signat_vert_comp(g):
    """
    """
    signat_vertices={}
    for vid in g.vertices():
        s=[]
        for vid2 in g.vertices():
            if vid2 in g.out_neighbors(vid):
                for eid in g.edges():
                    if g._edges[eid]==(vid,vid2):
                        s.append(g._edge_property['weight'][eid])
            else:
                s.append(0)
        signat_vertices[vid]=s                             
    return signat_vertices

def signat_edg_comp(g, signat_vertices):
    '''    
    '''
    signat_edges={}
    for eid in g.edges():
        vid1=g._edges[eid][0]
        vid2=g._edges[eid][1]
        s=signat_vertices[vid1][:]       
        s.append(s[vid2-1])
        s[vid2-1]=0
        signat_edges[eid]=s
    return signat_edges

def quasi_isomorphic_edges(eid1, eid2, signat_edges):
    '''
    '''
    if signat_edges[eid1]==signat_edges[eid2]:
        return True
    else:
        return False

def descendents(vid, g):
    '''
    '''
    d=[vid]
    for v in d:
        for v2 in g.out_neighbors(v):
            if v2 not in d:
                d.append(v2)
    d=d[1:]
    return d

def parents(vid,E):
    """D=edge.items()
    par=[]
    for i in range(len(D)):
        if vid in D[i][1]:
            par.append(D[i][0])"""
    par=[]
    K=E.keys()
    for i in range(len(K)):
        if K[i][1]==vid:
            par.append(K[i][0])
    return par
                         
def find_eid(vid1, vid2, g):
    '''
    '''
    for eid in g.edges():
        if g._edges[eid]==(vid1,vid2):
            break
    return eid

def Quasi_isom_desc_vertices(vid, g, signat_edges):
    '''
    '''
    QIDV=[]
    for vid2 in descendents(vid, g):
        for vid1 in g.out_neighbors(vid):
            eid=find_eid(vid, vid1, g)
            for vid21 in g.out_neighbors(vid2):
                eid2=find_eid(vid2, vid21, g)
                if quasi_isomorphic_edges(eid, eid2, signat_edges)==True:
                    if eid not in QIDV:
                        QIDV.append(eid)
                    QIDV.append(eid2)
                    break
    return QIDV

def MQPP(vid, g, edge, signat_edges):
    '''
    '''
    M=[]
    stop=False
    QIP=[]
    for eid in Quasi_isom_desc_vertices(vid, g, signat_edges):
        QIP.append([eid])
    while (len(QIP)>1) and stop==False:
        if len(QIP)==0:
            QPP=QIP[:]
        else:
            """construire QPP et tester si QPP est valable et donc arret"""
            QPP=[]
            QPP.append(QIP[0][:])
            Q=QIP[1:len(QIP)]
            a=True
            while a==True:
                elt=QPP[len(QPP)-1]
                eid=elt[len(elt)-1]
                i=0
                trouv=False
                while trouv==False and i<len(Q):
                    eid2=Q[i][0]
                    if g._edges[eid][1]==g._edges[eid2][0]:
                        QPP.append(Q[i])
                        Q=Q[0:i]+Q[i+1:]
                        trouv=True
                        a=True
                    else:
                        i=i+1
                if trouv==False:
                    a=False                        
            if len(QPP)>1:
                M.append(QPP)
                l=len(QPP)-1
                vidf=g._edges[QPP[l][len(QPP[l])-1]][1]
                D1=[]
                for eid in Quasi_isom_desc_vertices(vidf, g, signat_edges):
                    vid1=g._edges[eid][0]
                    vid2=g._edges[eid][1]
                    D1.append((vid1,vid2))
                find=False
                i=0
                while (not find) and i<len(Q):
                    lQ=len(Q[i])-1
                    vidf2=g._edges[Q[i][lQ]][1]
                    if vidf2 in edge.keys():
                        D2=edge[vidf2]
                        cmpt=0
                        for j in range(len(D2)):
                            if (vidf2,D2[j]) in D1:
                                cmpt+=1
                        if cmpt>0:
                            find=True
                    i+=1
                    
                if (len(QPP)==len(QIP)) or not find:
                    stop=True
        """ ralonger les paths de QIP"""
        if stop==False:
            vid2=g._edges[QIP[0][len(QIP[0])-1]][1]
            i=0
            while i<len(QIP):
                eid=QIP[i][len(QIP[i])-1]
                vid=g._edges[eid][1]
                add=False
                for eid2 in Quasi_isom_desc_vertices(vid2, g, signat_edges):
                    vid3=g._edges[eid2][0]
                    if vid==vid3:
                        QIP[i].append(eid2)
                        add=True
                        break
                if add==False:
                    QIP=QIP[0:i]+QIP[i+1:]
                else:
                    i=i+1

    return M

def paths_inclusion(p1, p2, g):
    '''
    '''
    P=[]
    for i in range(len(p1)):
        P=P+p1[i]
    Q=[]
    for i in range(len(p2)):
        Q=Q+p2[i]
    inclusion=False
    if len(P)<len(Q):
        cmp=0
        for i in range(len(P)):
            if P[i] in Q:
                cmp=cmp+1
        if cmp==len(P):
                inclusion=True
    return inclusion

def paths_perfect_inclusion(p1, p2):
    '''
    '''
    P=[]
    for i in range(len(p1)):
        P=P+p1[i]        
    perfect_inclusion=False
    if len(P)<len(p2[0]):
        for i in range(len(p2)):
            if P<=p2[i]:
                perfect_inclusion=True
                break
    return perfect_inclusion

def MQPPRT(g, edge, signat_edges):
    '''
    '''
    MRT=[]
    for vid in g.vertices():
        MM=MQPP(vid, g, edge, signat_edges)
        for j in range(len(MM)):
            M=MM[j]
            if len(M)>0:            
                ok=True
                for i in range(len(MRT)):
                    if paths_inclusion(M, MRT[i], g) and not(paths_perfect_inclusion(M, MRT[i])):        
                        ok=False
                if ok:
                    MRT.append(M)
                    i=0
                    while i<len(MRT)-1:
                        if paths_inclusion(MRT[i], M, g) and not(paths_perfect_inclusion(MRT[i], M)):
                            MRT.remove(MRT[i])
                        else:
                            i=i+1
    return MRT
        
def paths_intersection(p1, p2, g, list_node):
    '''
    '''
    m=[]
    for j in range(len(p1)):
        for k in range(len(p1[j])):
            eid=p1[j][k]
            vid1=g._edges[eid][0]
            vid2=g._edges[eid][1]
            if list_node[vid1-1] not in m:
                m.append(list_node[vid1-1])
            if list_node[vid2-1] not in m:
                m.append(list_node[vid2-1])

    n=[]
    for j in range(len(p2)):
        for k in range(len(p2[j])):
            eid=p2[j][k]
            vid1=g._edges[eid][0]
            vid2=g._edges[eid][1]
            if list_node[vid1-1] not in n:
                n.append(list_node[vid1-1])
            if list_node[vid2-1] not in n:
                n.append(list_node[vid2-1])

    intersect=False
    for i in range(len(m)):
        if m[i] in n:
            intersect=True
            break
    return intersect

def MQPPRT_classes(g, paths, list_node):
    '''
    '''
    MQPPRTINDEP=[]
    MQPPRTOVERL=[]
    indep=[]
    for i in range(len(paths)):
        indep.append(True)
    for i in range(len(paths)-1):
        k=0
        trouv=False
        while k<len(MQPPRTOVERL):
            if paths[i] not in MQPPRTOVERL[k]:
                k=k+1
            else:
                trouv=True
                break
        if k==0 or not trouv:
            MQPPRTOVERL.append([paths[i]])
        """else:
            MQPPRTOVERL[k].append(paths[i])"""
        
        j=i+1
        while j<len(paths):
            if paths_intersection(paths[i], paths[j], g, list_node)==True: 
                indep[i]=False
                indep[j]=False
                """if not(paths_perfect_inclusion(paths[i], paths[j]) or paths_perfect_inclusion(paths[j], paths[i])):"""
                if paths[j] not in MQPPRTOVERL[k]:
                    MQPPRTOVERL[k].append(paths[j])
            j+=1

    i=0
    while i<len(MQPPRTOVERL):
        if len(MQPPRTOVERL[i])==1:
            MQPPRTOVERL.remove(MQPPRTOVERL[i])
        else: i=i+1

    for i in range(len(indep)):
        if indep[i]==True:
            MQPPRTINDEP.append(paths[i])
    
    """ MQPPRTINCLUDED computation"""        
    p=[]
    for i in range(len(paths)):
        if paths[i] not in MQPPRTINDEP:
            p.append(paths[i])

    lenghts=[]
    for i in range(len(p)):
        l=0
        for j in range(len(p[i])):
            l=l+len(p[i][j])
        lenghts.append([l,i])
    lenghts.sort()
    sortedp=[]
    for i in range(len(lenghts)):
        sortedp.append(p[lenghts[i][1]])

    MQPPRTINCL=[]
    vue=[]
    for i in range(len(sortedp)):
        vue.append(False)
    for k in range(len(sortedp)-1):
        i=k
        if vue[i]==False:
            vue[i]=True
            inclset=[sortedp[i]]
            j=i+1
            while j<len(sortedp):
                if paths_perfect_inclusion(sortedp[i],sortedp[j]):
                    inclset.append(sortedp[j])
                    vue[i]=True
                    i=j
                j=j+1
            if len(inclset)>1:
                MQPPRTINCL.append(inclset)
                
    """supprimer les elt de overl qui sont en double dans incl"""

    k=0
    while k<len(MQPPRTOVERL):

        lenghts=[]
        for i in range(len(MQPPRTOVERL[k])):
            l=0
            for j in range(len(MQPPRTOVERL[k][i])):
                l=l+len(MQPPRTOVERL[k][i][j])
            lenghts.append([l,i])
        lenghts.sort()
        sortedp=[]
        for i in range(len(lenghts)):
            sortedp.append(MQPPRTOVERL[k][lenghts[i][1]])
        stop=False
        k1=0
        while k1<len(sortedp)-1 and not stop:
            if paths_perfect_inclusion(sortedp[k1],sortedp[len(sortedp)-1]):
                k1=k1+1
            else:
                stop=True
        if not stop:
            MQPPRTOVERL.remove(MQPPRTOVERL[k])
        else:
            k=k+1
            
    return MQPPRTINDEP, MQPPRTOVERL,MQPPRTINCL

def simple_path_wining(g,P):
    
    p1=[]
    for i in range(len(P)):
        p1=p1+P[i]
    N=len(p1)+1
    n1=len(P[0])+1
    w=N-n1
    return w
                
def path_wining(g, P, MQPPRTINCL):
    '''
    '''
    i=0
    trouv=False
    while i<len(MQPPRTINCL) and not trouv:
        j=0
        trouv=False
        while j<len(MQPPRTINCL[i]) and not trouv:
            if MQPPRTINCL[i][j]==P:
                trouv=True
            else:
                j+=1
        if not trouv:
            i+=1
    if not trouv:
        w=simple_path_wining(g,P)
    else:
        w=0
        k=0
        while k<=j:
           w=w+simple_path_wining(g,MQPPRTINCL[i][k])
           k+=1
        
    return w

def overlap_select(g, list_node, MQPPRTOVERL, MQPPRTINCL):
    '''
    '''
    OV=MQPPRTOVERL  
    select=[]
    for i in range(len(OV)):
        """construction du graphe des indépendances"""
        noeuds=[]
        for j in range(len(OV[i])):
            noeuds.append((j,path_wining(g,OV[i][j],MQPPRTINCL)))
        arcs=[]
        for j in range(len(OV[i])-1):
            k=j+1
            while k<len(OV[i]):
                if (paths_intersection(OV[i][j],OV[i][k],g,list_node)==False) or (paths_perfect_inclusion(OV[i][j],OV[i][k]) or paths_perfect_inclusion(OV[i][k],OV[i][j])):
                    arcs.append((j,k))
                k=k+1
        """rechercher la clique max de gain max dans ce graphe"""
        cliques=[]
        for j in range(len(noeuds)):
            cliques.append([noeuds[j]])
        step=1
        stop=False
        while step<len(noeuds) and stop==False:
            print 'CLIQUES',cliques
            stop=True
            for j in range(len(cliques)):
                sets=[]
                for k in range(len(noeuds)):
                    if noeuds[k] not in cliques[j]:
                        sets.append(noeuds[k])
                add=[]
                for k in range(len(sets)):
                    cpt=0
                    for k1 in range(len(cliques[j])):
                        if ((sets[k],cliques[j][k1]) in arcs) or ((cliques[j][k1],sets[k]) in arcs):
                            cpt=cpt+1
                    if cpt==len(cliques[j]):
                        add.append(sets[k])
                        stop=False
                cliques[j]=cliques[j]+add
            """supprimer les cliques en double"""
            CL=[]
            for j in range(len(cliques)):
                if cliques[j] not in CL:
                    CL.append(cliques[j])
            cliques=CL
        """selectionner la clique de gain max"""
        gain=[]
        for j in range(len(cliques)):
            cpt=0
            for k in range(len(cliques[j])):
                trouv=False
                k1=0
                while k1<len(cliques[j]) and not(trouv):
                    if k1<>k:
                        if paths_perfect_inclusion(OV[i][cliques[j][k][0]],OV[i][cliques[j][k1][0]]):
                            trouv=True
                    k1=k1+1
                if not trouv:
                    cpt=cpt+cliques[j][k][1]
            gain.append(cpt)
        maxc=gain[0]
        j=1
        ind=0
        while j<len(gain):
            if gain[j]>=maxc:
                maxc=gain[j]
                ind=j
            j=j+1
        sol=[]
        for j in range(len(cliques[ind])):
            sol.append(OV[i][cliques[ind][j][0]])
        select.append(sol)        
    return select

""" Gre derivation"""

def DAG_init(g):
    node =[]
    for vid in g.vertices():
        node.append(vid)
    edge={}
    for vid in g.vertices():
        child=[]
        for eid in g.edges():
            if g._edges[eid][0]==vid:
                child.append(g._edges[eid][1])
        if len(child)>0:
            edge[vid]=child
    val={}
    for eid in g.edges():
        val[g._edges[eid]]=g._edge_property['weight'][eid]    
    return node,edge,val
        
   
def quotient_graph(g,node,edge,val,MQPPRTINDEP,MQPPRTINCL):
    '''
    '''
    NBsupE=0
    V=node[:]
    E={}
    for i in range(len(val.keys())):
        E[val.keys()[i]]=val[val.keys()[i]]
    RE={}
    FE={}
    PE={}
    initial_returns=[]

    for i in range(len(MQPPRTINDEP)):
        gre=reduceindep(g,node,edge,val,NBsupE,V,E,RE,FE,PE,MQPPRTINDEP[i])
        V=gre[0]
        E=gre[1]
        RE=gre[2]
        FE=gre[3]
        PE=gre[4]
        if gre[5]<>[]:
            initial_returns=gre[5]
        NBsupE=gre[6]

    for i in range(len(MQPPRTINCL)):
        gre=reduceincl(g,node,edge,val,NBsupE,V,E,RE,FE,PE,MQPPRTINCL[i])
        V=gre[0]
        E=gre[1]
        RE=gre[2]
        FE=gre[3]
        PE=gre[4]
        if gre[5]<>[]:
            initial_returns=gre[5]
        NBsupE=gre[6]
 
    return V,E,RE,FE,PE,initial_returns,NBsupE

def reduceindep(g,node,edge,val,NBsupE,V,E,RE,FE,PE,path):
    # a simple path reduction
    """path transformation from edges lists to nodes lists"""
    P=[]
    for i in range(len(path)):
        m=[]
        for j in range(len(path[i])):
            eid=path[i][j]
            vid1=g._edges[eid][0]
            vid2=g._edges[eid][1]
            if node[vid1-1] not in m:
                m.append(node[vid1-1])
            if node[vid2-1] not in m:
                m.append(node[vid2-1])
        P.append(m)
    print 'path',P

    """create the return edge"""
    RE[(P[0][len(P[0])-2],P[0][0])]=val[(P[0][len(P[0])-2],P[0][len(P[0])-1])]
    
    """create the final edge"""
    FE[(P[0][len(P[0])-2],P[len(P)-1][len(P[len(P)-1])-1])]=val[(P[0][len(P[0])-2],P[0][len(P[0])-1])]
    if (P[0][len(P[0])-2],P[0][len(P[0])-1]) in E.keys():
        E.pop((P[0][len(P[0])-2],P[0][len(P[0])-1]))

    """supression des arcs du path P"""
    i=1
    while i<len(P):
        for j in range(len(P[i])-1):
            if (P[i][j],P[i][j+1]) in E.keys():
                E.pop((P[i][j],P[i][j+1]))
                NBsupE+=1
        i+=1

    """redirect the parents edges to the principal path of P"""
    pp=[]
    for i in range(len(P)):
        for j in range(len(P[i])):
            pp.append(P[i][j])
    for i in range(len(P)):
        alpha=len(P)-(i+1)
        for j in range(len(P[i])-1):
            par=parents(P[i][j],E)
            for k in range(len(par)):
                if par[k] not in pp:
                    PE[(par[k],P[0][j],val[(par[k],P[i][j])])]=[alpha]
                    if (par[k],P[i][j]) in E.keys():
                        E.pop((par[k],P[i][j]))
    initial_returns=[]
    if len(parents(P[0][0],E))==0 and len(parents(P[0][0],PE))==0 and len(parents(P[0][0],FE))==0:
        initial_returns.append(len(P)-1)

    """suppression des sommets inaccessibles de V"""
    if True:
        stop=False
        while not stop:
            stop=True
            KE=[]
            for i in range(len(E.keys())):
                KE.append(E.keys()[i][1])
            KRE=[]
            for i in range(len(RE.keys())):
                KRE.append(RE.keys()[i][1])
            KFE=[]
            for i in range(len(FE.keys())):
                KFE.append(FE.keys()[i][1])
            KPE=[]
            for i in range(len(PE.keys())):
                KPE.append(PE.keys()[i][1])
            i=0
            while i<len(V):
                if (V[i]<>g.root) and (V[i] not in KE) and (V[i] not in KRE) and (V[i] not in KFE) and (V[i] not in KPE):
                    stop=False
                    """supprimer les edges sortants de V[i]"""
                    j=0
                    while j<len(E.keys()):
                        if E.keys()[j][0]==V[i]:
                            E.pop(E.keys()[j])
                            NBsupE+=1
                        else:
                            j+=1
                    j=0
                    while j<len(PE.keys()):
                        if PE.keys()[j][0]==V[i]:
                            PE.pop(PE.keys()[j])
                            NBsupE+=1
                        else:
                            j+=1
                    V.remove(V[i])
                else:
                    i+=1
        
        K=PE.keys()
        for i in range(len(K)):
            if (K[i][0] not in V) or (K[i][1] not in V):
                PE.pop(K[i])
    
    return V,E,RE,FE,PE,initial_returns,NBsupE

def reduceincl(g,node,edge,val,NBsupE,V,E,RE,FE,PE,seq):

    a=len(seq)-1
    initial_returns=[]
    while a>=0:
        """path transformation from edges lists to nodes lists"""
        P=[]
        path=seq[a]
        for i in range(len(path)):
            m=[]
            for j in range(len(path[i])):
                eid=path[i][j]
                vid1=g._edges[eid][0]
                vid2=g._edges[eid][1]
                if node[vid1-1] not in m:
                    m.append(node[vid1-1])
                if node[vid2-1] not in m:
                    m.append(node[vid2-1])
            P.append(m)
        print 'path',P
        
        """create the return edge"""
        RE[(P[0][len(P[0])-2],P[0][0])]=val[(P[0][len(P[0])-2],P[0][len(P[0])-1])]

        """create the final edge"""
        FE[(P[0][len(P[0])-2],P[len(P)-1][len(P[len(P)-1])-1])]=val[(P[0][len(P[0])-2],P[0][len(P[0])-1])]
        if (P[0][len(P[0])-2],P[0][len(P[0])-1]) in E.keys():
            E.pop((P[0][len(P[0])-2],P[0][len(P[0])-1]))

        """supression des arcs du path P"""
        i=1
        while i<len(P):
            for j in range(len(P[i])-1):
                if (P[i][j],P[i][j+1]) in E.keys():
                    E.pop((P[i][j],P[i][j+1]))
                    NBsupE+=1
            i+=1
            
        """redirect the parents edges"""
        pp=[]
        for i in range(len(P)):
            for j in range(len(P[i])):
                pp.append(P[i][j])
        for i in range(len(P)):
            alpha=len(P)-(i+1)
            for j in range(len(P[i])-1):
                if a==len(seq)-1:
                    par=parents(P[i][j],E)
                    for k in range(len(par)):
                        if par[k] not in pp:
                            PE[(par[k],P[0][j],val[(par[k],P[i][j])])]=[alpha]
                            if (par[k],P[i][j]) in E.keys():
                                E.pop((par[k],P[i][j]))
                else:
                    par=PE.keys()
                    for k in range(len(par)):
                        if par[k][1]==P[i][j]:
                            PE[(par[k][0],P[0][j],par[k][2])]=PE[(par[k][0],P[i][j],par[k][2])]
                            PE[(par[k][0],P[0][j],par[k][2])].append(alpha)
                            if P[i][j]<>P[0][j]:
                                PE.pop((par[k][0],P[i][j],par[k][2]))
        a=a-1
        if len(parents(P[0][0],E))==0 and len(parents(P[0][0],PE))==0 and len(parents(P[0][0],FE))==0: 
            initial_returns.append(len(P)-1)

        """suppression des sommets superflus"""
        stop=False
        while not stop:
            stop=True
            KE=[]
            for i in range(len(E.keys())):
                KE.append(E.keys()[i][1])
            KRE=[]
            for i in range(len(RE.keys())):
                KRE.append(RE.keys()[i][1])
            KFE=[]
            for i in range(len(FE.keys())):
                KFE.append(FE.keys()[i][1])
            KPE=[]
            for i in range(len(PE.keys())):
                KPE.append(PE.keys()[i][1])
            i=0
            while i<len(V):
                if (V[i]<>g.root) and (V[i] not in KE) and (V[i] not in KRE) and (V[i] not in KFE) and (V[i] not in KPE):
                    stop=False
                    """supprimer les edges sortants de V[i]"""
                    j=0
                    while j<len(E.keys()):
                        if E.keys()[j][0]==V[i]:
                            E.pop(E.keys()[j])
                            NBsupE+=1
                        else:
                            j+=1
                    j=0
                    while j<len(PE.keys()):
                        if PE.keys()[j][0]==V[i]:
                            PE.pop(PE.keys()[j])
                            NBsupE+=1
                        else:
                            j+=1
                    V.remove(V[i])
                else:
                    i+=1

        K=PE.keys()
        for i in range(len(K)):
            if (K[i][0] not in V) or (K[i][1] not in V):
                PE.pop(K[i])

    return V,E,RE,FE,PE,initial_returns,NBsupE

if(__name__ == "__main__"):
    
    """mtg_filename = './example_MTG.mtg'
    mtg_filename = './fractal00.mtg'

    g = MTG(mtg_filename)

    print 'Create the graph'
    tree = mtg2graph(1, 2)

    print 'Reduction'
    dag_exact = tree_reduction(tree)
    print is_linear(dag_exact)
    
    dag_linearization(dag_exact)
    
    tree_rec = tree_reconstruction(tree)
    
    
    print_graph(dag_exact)
    
    print 'nb-signature',dag_exact.nb_signature
    print 'max_signature',dag_exact.max_signature_edge
    print 'sum_signature_edge',dag_exact.sum_signature_edge
    print 'nodes'
    for i in range(len(dag_exact.nodes)):
        print dag_exact.nodes[i], parent_list(dag_exact.nodes[i])
    print 'edges',len(dag_exact.edges)
    for i in range(len(dag_exact.edges)):
        print dag_exact.edges[i].label, dag_exact.edges[i].start.signature, dag_exact.edges[i].end.signature


    print dag_exact.set_nb_signature(nb_nodes(dag_exact))
    
    print is_linear(dag_exact)
    G=mapping_cgraph_test1(dag_exact)"""
    
    G=test1()
    #print 'taille',taille_arbre(G,1)
    List1=DAG_init(G)
    list_node=List1[0]
    edge=List1[1]
    val=List1[2]
    S=signat_vert_comp(G)
    E=signat_edg_comp(G, S)
    Paths=MQPPRT(G, edge, E)
    print 'paths', Paths
    MQPPRTINDEP=MQPPRT_classes(G, Paths, list_node)[0]
    MQPPRTOVERL=MQPPRT_classes(G, Paths, list_node)[1]
    MQPPRTINCL=MQPPRT_classes(G, Paths, list_node)[2]
    print 'MQPPRTINDEP',MQPPRTINDEP
    print 'MQPPOVERL',MQPPRTOVERL
    print 'MQPPINCL',MQPPRTINCL
    selection=overlap_select(G,list_node,MQPPRTOVERL,MQPPRTINCL)
    print 'selection', selection
    for i in range(len(MQPPRTOVERL)):
        for j in range(len(MQPPRTOVERL[i])):
            p=MQPPRTOVERL[i][j]
            if p in selection[i]:
                k=0
                trouv=False
                while (trouv==False) and (k<len(MQPPRTINCL)):
                    if p in MQPPRTINCL[k]:
                        trouv=True
                    k=k+1
                if trouv==False:
                    MQPPRTINDEP.append(p)
            else:
                for k in range(len(MQPPRTINCL)):
                    if p in MQPPRTINCL[k]:
                        MQPPRTINCL[k].remove(p)
    i=0
    while i<len(MQPPRTINCL):
        if len(MQPPRTINCL[i])==1:
            MQPPRTINDEP.append(MQPPRTINCL[i][0])
            MQPPRTINCL.remove(MQPPRTINCL[i])
        else:
            i+=1

    INC=[]
    for i in range(len(MQPPRTINCL)):
        vrai=True
        j=0
        while vrai and j<len(MQPPRTINCL)-1:
            P=[]
            for k in range(len(MQPPRTINCL[i][j])):
                P=P+MQPPRTINCL[i][j][k]
            Q=MQPPRTINCL[i][j+1][0]
            inclusion=False
            if len(P)<len(Q):
                cmp=0
                for k in range(len(P)):
                    if P[k] in Q:
                        cmp=cmp+1
                if cmp==len(P):
                        inclusion=True                
            if not inclusion:
                vrai=False
            else:
                j+=1
        if vrai:
            INC.append(MQPPRTINCL[i])
    MQPPRTINCL=INC
    print 'MQPPRTINDEP',MQPPRTINDEP
    print 'MQPPOVERL',MQPPRTOVERL
    print 'MQPPINCL',MQPPRTINCL
    List=quotient_graph(G,list_node,edge,val,MQPPRTINDEP,MQPPRTINCL)
    V=List[0]
    E=List[1]
    RE=List[2]
    FE=List[3]
    PE=List[4]
    initial_returns=List[5]
    NBsupE=List[6]
    print 'resulting graph'
    print 'V',V
    print 'E',E
    print 'RE',RE
    print 'FE',FE
    print 'PE',PE
    print 'initial returns',initial_returns
            
    TN=len(list_node)
    Tedge=len(val.keys())
    initial=TN+Tedge
    print initial
    TV=len(V)
    TE=len(E.keys())
    TRE=len(RE.keys())
    TFE=len(FE.keys())
    TPE=len(PE.keys())
    final=TV+TE+TRE+TFE+TPE
    print final
    if (initial>0):
        gain=initial-final
        print 'the gain',gain
        icf=1-(final.__float__()/initial.__float__())
        print 'the improved compression factor=',icf
        #print nRe.__float__(), nR.__float__()
