# -*- coding: cp1252 -*-

from amlPy import *

from vplants.self_similarity.graph import *
from vplants.self_similarity.cst import *
from vplants.self_similarity.mtg2graph import *

from vplants.tree_reduction.graph import *
from vplants.tree_reduction.functions import *
from openalea.container.graph import *
import string



"""   rec directe
def test1():
    edges=[(1,2),(2,3),(3,4),(3,5),(4,6),(4,5),(5,6)]
    edges_weight={0:1,1:2,2:1,3:1,4:1,5:1,6:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 6
    assert g.nb_edges() == 7
    return g"""

""" rec indirecte
def test1():
    edges=[(1,2),(1,7),(2,3),(2,6),(3,4),(3,7),(4,6),(4,5),(5,7),(6,8),(7,8)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:3,10:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 8
    assert g.nb_edges() == 11
    return g"""

""" perfectly included paths
def test1():
    edges=[(1,2),(1,7),(2,3),(2,7),(3,4),(3,8),(4,7),(4,5),(5,6),(5,7),(6,8),(6,9),(8,9),(7,9)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:2,13:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 9
    assert g.nb_edges() == 14
    return g"""

""" perfectly included paths
def test1():
    edges=[(1,2),(1,7),(2,3),(2,7),(3,4),(3,8),(4,7),(4,5),(5,6),(5,7),(6,8),(6,9),(8,9),(7,9)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:2,13:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 9
    assert g.nb_edges() == 14
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

"""   T0,4
def test1():
    edges=[(1,2),(2,3),(2,4),(3,5),(3,6),(4,7),(4,15),(5,8),(5,17),(6,9),(6,17),(7,10),(7,15),(8,11),(8,13),(9,12),(9,15),(10,15),(11,14),(11,17),(12,17),(13,16),(14,18),(14,22),(15,19),(15,22),(16,21),(17,19),(18,20),(18,22),(19,21),(19,22),(20,22),(21,22),(22,23)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1,15:1,16:1,17:2,18:1,19:1,20:2,21:1,22:1,23:1,24:1,25:1,26:1,27:1,28:1,29:1,30:1,31:1,32:2,33:2,34:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 23
    assert g.nb_edges() == 35
    return g"""

"""   T0,8
def test1():
    edges=[(1,2),(2,3),(3,4),(4,5),(4,6),(5,7),(5,9),(6,8),(7,9),(8,10),(8,13),(9,11),(9,13),(10,13),(11,12),(12,13),(13,14)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:2,9:1,10:1,11:1,12:1,13:2,14:1,15:2,16:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 14
    assert g.nb_edges() == 17
    return g"""

"""    T1
def test1():
    edges=[(1,2),(1,4),(2,3),(2,4),(3,4),(4,5),(4,7),(5,6),(5,7),(6,7),(7,8),(7,10),(8,9),(8,10),(9,10),(10,11),(11,12),(12,13),(12,14),(13,14),(14,15),(14,16),(15,16),(16,17)]
    edges_weight={0:1,1:1,2:1,3:1,4:2,5:1,6:1,7:1,8:1,9:2,10:1,11:1,12:1,13:1,14:2,15:1,16:2,17:1,18:1,19:2,20:1,21:1,22:2,23:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 17
    assert g.nb_edges() == 24
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
    edges=[(1,2),(2,3),(3,4),(4,5),(5,6),(6,7)]
    edges_weight={0:3,1:3,2:3,3:3,4:3,5:4}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 7
    assert g.nb_edges() == 6
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
    return g
"""

"""    Arbre1 séquence 
def test1():
    edges=[(1,2),(2,3),(3,4)]
    edges_weight={0:1,1:1,2:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 4
    assert g.nb_edges() == 3
    return g"""

"""    Arbre2 séquence 
def test1():
    edges=[(1,2),(2,3),(2,4),(3,4),(4,5)]
    edges_weight={0:1,1:1,2:1,3:1,4:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 5
    assert g.nb_edges() == 5
    return g
"""
"""    Arbre3 séquence 
def test1():
    edges=[(1,2),(2,3),(2,5),(3,4),(3,5),(4,5),(5,6)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 6
    assert g.nb_edges() == 7
    return g
"""
"""    Arbre4 séquence 
def test1():
    edges=[(1,2),(2,3),(2,6),(3,4),(3,6),(4,5),(4,6),(5,6),(6,7)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 7
    assert g.nb_edges() == 9
    return g
"""
"""    Arbre augmenté séquence """
def test1():
    edges=[(1,2),(1,5),(1,9),(1,14),(2,3),(3,4),(4,20),(5,6),(6,7),(6,8),(7,8),(8,20),(9,10),(10,11),(10,13),(11,12),(11,13),(12,13),(13,20),(14,15),(15,16),(15,19),(16,17),(16,19),(17,18),(17,19),(18,19),(19,20)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:3,7:1,8:1,9:1,10:1,11:3,12:1,13:1,14:1,15:1,16:1,17:1,18:3,19:1,20:1,21:1,22:1,23:1,24:1,25:1,26:1,27:3}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 20
    assert g.nb_edges() == 28
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

def MQPP(vid, g, signat_edges):
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
            Q=QIP[:]
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
                M=M+QPP
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

def MQPPRT(g, signat_edges):
    '''
    '''
    MRT=[]
    for vid in g.vertices():
        M=MQPP(vid, g, signat_edges)
        ok=True
        for i in range(len(MRT)):
            if paths_inclusion(M, MRT[i], g):
                ok=False
        if ok:
            MRT.append(M)
            i=0
            while i<len(MRT)-1:
                if paths_inclusion(MRT[i], M, g):
                    MRT.remove(MRT[i])
                else:
                    i=i+1         
    return MRT
        
def paths_intersection(p1, p2, g):
    '''
    '''
    P=[]
    for i in range(len(p1)):
        P=P+p1[i]
    Q=[]
    for i in range(len(p2)):
        Q=Q+p2[i]
    intersect=False
    for i in range(len(P)):
        if P[i] in Q:
            intersect=True
            break
    return intersect

def MQPPRT_classes(g, paths):
    '''
    '''
    MQPPRTINDEP=[]
    MQPPRTOVERL=[]
    indep=[]
    for i in range(len(paths)):
        indep.append(True)
    for i in range(len(paths)-1):
        j=i+1
        while j<len(paths):
            if paths_intersection(paths[i], paths[j], g)==True:
                MQPPRTOVERL.append((paths[i], paths[j]))
                indep[i]=False
                indep[j]=False
            j+=1
    for i in range(len(indep)):
        if indep[i]==True:
            MQPPRTINDEP.append(paths[i])

    return MQPPRTINDEP, MQPPRTOVERL
                
def path_wining(g, P):
    '''
    '''
    p1=[]
    for i in range(len(P)):
        p1=p1+P[i]
    N=len(p1)+1
    n1=len(P[0])+1
    w=N-n1
    return w

def overlap_select(g, MQPPRTOVERL):
    '''
    '''
    select=[]
    for i in range(len(MQPPRTOVERL)):
        w1=path_wining(g, MQPPRTOVERL[i][0])
        w2=path_wining(g, MQPPRTOVERL[i][1])
        if w1<w2:
            if MQPPRTOVERL[i][1] not in select:
                select.append(MQPPRTOVERL[i][1])
        elif MQPPRTOVERL[i][0] not in select:
            select.append(MQPPRTOVERL[i][0])

    return select

"""GRT creation and GRTLre derivation"""

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
        
def G_construction(node_list,edge_list,val_edges):
    """ a partir du DAG construire le graphe G équivalent, G est représenté par:
    list_node=['a','b',...] l'ensemble des symboles associés aux nœuds de G
    G_node={'a':['chaine de a'],....} chaque nœuds dispose de la liste des chaines qui lui sont attribués, cette liste contient au max 2 éléments (en cas de récursivité)
    list_node et G_node représentent les informations concernant les nœuds
    list_edge=[[list_node.index('a'),[list_node.index(des fils de 'a'),...]],...] pour chaque nœud de list_node donner la liste de ces nœuds fils
    G_edge={('a','fils de a'):chaine associé a l'arc en question,... } donne pour chaque arc la chaine correspondante
    list_edge et G_edge représentent les informations concernant les arcs.
    Une structure supplémentaire est nécessaire pour cette application,
    c'est nbr_edge elle contient pour chaque valeur d'arc (element de G_edge.vals) son nombre d'apparitions dans l'arbre
    nbr_edge={chaine1:nbr,....} """
    list_node=[]
    list_edge=[]
    G_node={}
    G_edge={}
    equiv={}
    parameters={}
    parameters1={}
    parameters2={}
    V=[]
    for i in range(len(node_list)):
        symbol=creer_symbol(V)
        V=V+[symbol]
        equiv[node_list[i]]=symbol
    for i in range(len(node_list)):
        sb=''
        if i==len(node_list)-1:
            sb='l'
        else:
            for j in range(len(edge_list[node_list[i]])):
                fils=edge_list[node_list[i]][j]
                nbr=val_edges[(node_list[i],fils)]
                for k in range(nbr):
                    sb=sb+'[F'+equiv[fils]+']'                                
        G_node[equiv[node_list[i]]]=[sb]
        list_node=list_node+[equiv[node_list[i]]]
    for i in range(len(node_list)-1):      
        list_edge=list_edge+[[i,[]]]
        for j in range(len(edge_list[node_list[i]])):
            k1=0
            while (node_list[k1]!=edge_list[node_list[i]][j]):
              k1=k1+1
            list_edge[i][1]=list_edge[i][1]+[k1]
            G_edge[(equiv[node_list[i]],equiv[edge_list[node_list[i]][j]])]=string.replace(G_node[equiv[node_list[i]]][0],equiv[edge_list[node_list[i]][j]],'W')
    nbr_edge=nbr_arc_construction(G_edge)
    for k in range(len(G_edge.keys())):
      for j in range(len(nbr_edge.keys())):
        if equal_string(G_edge[G_edge.keys()[k]],nbr_edge.keys()[j]):
          G_edge[G_edge.keys()[k]]=nbr_edge.keys()[j]    
    return list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2
    

def quotient_graph(g,list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2, setP):
    '''
    '''
    setPg=[]
    for i in range(len(setP)):
        P=setP[i]
        Pg=[]
        for j in range(len(P)):
            m=[]
            for k in range(len(P[j])):
                eid=P[j][k]
                vid1=g._edges[eid][0]
                vid2=g._edges[eid][1]
                if list_node[vid1-1] not in m:
                    m.append(list_node[vid1-1])
                if list_node[vid2-1] not in m:
                    m.append(list_node[vid2-1])
            Pg.append(m)
        setPg.append(Pg)

    for i in range(len(setPg)):
        """calcul du graphe quotient"""
        if (len(setPg[i])>0):
            L=treat_long_path(setPg[i],list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2)
            list_node=L[0]
            list_edge=L[1]
            G_node=L[2]
            G_edge=L[3]
            nbr_edge=L[4]
            parameters=L[5]
            parameters1=L[6]

    return list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2

def treat_long_path(glist,list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2):
    """elle reçoit en entrée un chemin de type
    glist=[[x1,x2,...xk],[xk,xk+1,...,x2k-1],[x2k-1,...,x3k-2],...,[xend,...x2end-1]]
    et effectue la réduction adéquate de G"""
    """modifier le noeud associé a xk par la valeur [x1,x2end-1]"""
    print glist
    k=len(glist[0])-2
    string1=G_node[glist[0][k]][0]
    string2=string.replace(string1, glist[0][k+1], glist[0][0])
    string3=string.replace(string1, glist[0][k+1], glist[len(glist)-1][len(glist[len(glist)-1])-1])
    G_node[glist[0][k]]=[string2,string3]
    parameters[glist[0][0]]=[glist[0][0],len(glist)-1]
    parameters2[glist[0][k]]=glist[0][0]
    l2=sup_edge(list_node.index(glist[0][k]),list_node.index(glist[0][k+1]),list_edge,list_node,G_edge,nbr_edge)
    list_edge=l2[0]
    G_edge=l2[1]
    nbr_edge=l2[2]

   
    """créer deux arcs entre xk-1 et x2end-1 et x1 ayant la valeur 'alpha'"""
    G_edge[(glist[0][k],glist[len(glist)-1][len(glist[len(glist)-1])-1])]='alpha'
    G_edge[(glist[0][k],glist[0][0])]='alpha'
    k1=ind_find(list_node.index(glist[0][k]),list_edge)
    if list_node.index(glist[len(glist)-1][len(glist[len(glist)-1])-1]) not in list_edge[k1][1]:
        list_edge[k1][1].append(list_node.index(glist[len(glist)-1][len(glist[len(glist)-1])-1]))
    if list_node.index(glist[0][0]) not in list_edge[k1][1]:
        list_edge[k1][1].append(list_node.index(glist[0][0]))

    """Parcourir les nœuds xi de xend jusqu'à x2k-1,
    si xi possède un seul nœud parent qui est xi-1 alors supprimer le lien
    de xi avec xi+1, sinon arrêter le parcours et effectuer le traitement (a) """
    i=len(glist)-1
    while (i>0):
        for t in range(len(glist[i])-1):
            if (nbr_parent(list_node.index(glist[i][t]),list_edge)>1):            
                """relier chaque parent xj de xi avec x1"""
                l_parent=list_parent(list_node.index(glist[i][t]),list_edge)
                k=0
                while (k<len(l_parent)):                
                    if l_parent[k]!=list_node.index(glist[i-1][len(glist[i-1])-2]):
                        """relier xj à x1 et supprimer le lien avec xi"""
                        k1=ind_find(l_parent[k],list_edge)    
                        list_edge[k1][1].remove(list_node.index(glist[i][t]))
                        if (nbr_edge[G_edge[(list_node[l_parent[k]],glist[i][t])]]>=1):
                            nbr_edge[G_edge[(list_node[l_parent[k]],glist[i][t])]]-=1
                        Y=G_edge[(list_node[l_parent[k]],glist[i][t])]                    
                        G_edge.pop((list_node[l_parent[k]],glist[i][t]))                                                     
                        list_edge[k1][1].append(list_node.index(glist[0][t]))
                        G_edge[(list_node[l_parent[k]],glist[0][t])]=Y                    
                        parameters[glist[i][t]]=[glist[0][t],len(glist)-i-1]
                    k+=1
        """supprimer le lien de chaque x dans [xi...xi+(k-1)-1] avec tous ces fils"""
        j=len(glist[i])-2
        while j>=0:
            k=ind_find(list_node.index(glist[i][j]),list_edge)            
            k1=0
            while (k1<len(list_edge[k][1])):
                l2=sup_edge(list_node.index(glist[i][j]),list_edge[k][1][k1],list_edge,list_node,G_edge,nbr_edge)
                list_edge=l2[0]
                G_edge=l2[1]
                nbr_edge=l2[2]
                k1+=1
            j-=1
        l2=sup_edge(list_node.index(glist[1][0]),list_node.index(glist[1][1]),list_edge,list_node,G_edge,nbr_edge)
        list_edge=l2[0]
        G_edge=l2[1]
        nbr_edge=l2[2]
        i-=1
    """parametrer les régles de glist[0]"""
    for i in range(len(glist[0])-2):
        parameters1[glist[0][i]]=glist[0][i+1]
        
    """Après la modification de G il faut supprimer les nœuds superflus (inaccessibles)"""
    l1=ajust_sup(list_edge,G_edge,list_node,G_node,nbr_edge)
    list_edge=l1[0]
    G_edge=l1[1]
    list_node=l1[2]
    G_node=l1[3]
    nbr_edge=l1[4]   
    return list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2

def L_system(G_node,list_node):
    """affichage du PD0L_système a partir du graphe G"""
    print 'Axiom=',list_node[0]
    print 'rules:'
    for i in range(len(list_node)):
        for j in range(len(G_node[list_node[i]])):
            print list_node[i],'->',G_node[list_node[i]][j]

def Parametric_L_system(G_node,list_node,parameters,parameters1,parameters2):  
    for k1 in range(len(G_node.keys())):
        if len(G_node[G_node.keys()[k1]][0])>1:
            for k2 in range(len(G_node[G_node.keys()[k1]])):
                for k in range(len(parameters.keys())):                
                    if parameters.keys()[k]!= G_node.keys()[k1] and parameters.keys()[k]==parameters[parameters.keys()[k]][0] and ((G_node.keys()[k1] not in parameters1.keys()) or ((G_node.keys()[k1] in parameters1.keys()) and (parameters.keys()[k]!=parameters1[G_node.keys()[k1]]))) and ((G_node.keys()[k1] not in parameters2.keys()) or ((G_node.keys()[k1] in parameters2.keys()) and (parameters.keys()[k]!=parameters2[G_node.keys()[k1]]))):
                        G_node[G_node.keys()[k1]][k2]=string.replace(G_node[G_node.keys()[k1]][k2],parameters.keys()[k],parameters[parameters.keys()[k]][0]+'('+str(parameters[parameters.keys()[k]][1])+')')
                for k in range(len(parameters.keys())):
                    if parameters.keys()[k]!= G_node.keys()[k1] and parameters.keys()[k]!=parameters[parameters.keys()[k]][0] and ((G_node.keys()[k1] not in parameters1.keys()) or ((G_node.keys()[k1] in parameters1.keys()) and (parameters.keys()[k]!=parameters1[G_node.keys()[k1]]))) and ((G_node.keys()[k1] not in parameters2.keys()) or ((G_node.keys()[k1] in parameters2.keys()) and (parameters.keys()[k]!=parameters2[G_node.keys()[k1]]))):                     
                        G_node[G_node.keys()[k1]][k2]=string.replace(G_node[G_node.keys()[k1]][k2],parameters.keys()[k],parameters[parameters.keys()[k]][0]+'('+str(parameters[parameters.keys()[k]][1])+')')
                    
    """affichage du L_système a partir du graphe G"""
    symbol=creer_symbol(G_node)
    if list_node[0] in parameters.keys():
        Axiom=parameters[list_node[0]][0]+'('+str(parameters[list_node[0]][1])+')'
    else:
        Axiom=list_node[0]
    print 'Axiom=',Axiom
    print 'rules:'
    for k1 in range(len(G_node.keys())):
        for k2 in range(len(G_node[G_node.keys()[k1]])):
            G_node[G_node.keys()[k1]][k2]=string.replace(G_node[G_node.keys()[k1]][k2],list_node[len(list_node)-1],'l')
        
    for i in range(len(list_node)-1):
        if (len(G_node[list_node[i]])==1):
            if list_node[i] in parameters1.keys():
                print list_node[i],'(',symbol,') ->',string.replace(G_node[list_node[i]][0],parameters1[list_node[i]],parameters1[list_node[i]]+'('+symbol+')')
            else:
                print list_node[i],'->',G_node[list_node[i]][0]
        else:
            if list_node[i] in parameters2.keys():
                print list_node[i],'(',symbol,')',':',symbol,'>0 ->',string.replace(G_node[list_node[i]][0],parameters2[list_node[i]],parameters2[list_node[i]]+'('+symbol+'-1)')
                print list_node[i],'(',symbol,')',':',symbol,'=0 ->',G_node[list_node[i]][1]

if(__name__ == "__main__"):
    """
    mtg_filename = './example_MTG.mtg'

    g = MTG(mtg_filename)

    print 'Create the graph'
    tree = mtg2graph(1, 2)

    print 'Reduction'
    dag_exact = tree_reduction(tree)
    tree_rec = tree_reconstruction(tree)"""
    G=test1()
    S=signat_vert_comp(G)
    E=signat_edg_comp(G, S)
    Paths=MQPPRT(G, E)
    MQPPRTINDEP=MQPPRT_classes(G, Paths)[0]
    MQPPRTOVERL=MQPPRT_classes(G, Paths)[1]
    selection=overlap_select(G,MQPPRTOVERL)
    set_paths=MQPPRTINDEP+selection
    List1=DAG_init(G)
    List=G_construction(List1[0],List1[1],List1[2])
    list_node=List[0]
    list_edge=List[1]
    G_node=List[2]
    G_edge=List[3]
    nbr_edge=List[4]
    parameters=List[5]
    parameters1=List[6]
    parameters2=List[7]
    nR=len(list_node)
    print 'the corresponding PD0L-system:'
    L_system(G_node,list_node)
    List=quotient_graph(G,list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1,parameters2,set_paths)
    list_node=List[0]
    list_edge=List[1]
    G_node=List[2]
    G_edge=List[3]
    nbr_edge=List[4]
    parameters=List[5]
    parameters1=List[6]
    parameters2=List[7]
    """affichage du L_système extrait de G"""
    print 'the corresponding parametric PD0L-system:'
    Parametric_L_system(G_node,list_node,parameters,parameters1,parameters2)
    nRe=len(list_node)
    icf=1-(nRe.__float__()/nR.__float__())
    print 'the improved compression factor=',icf

   
