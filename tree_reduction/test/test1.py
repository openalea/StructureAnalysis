# -*- coding: cp1252 -*-

"""from amlPy import *

from vplants.self_similarity.graph import *
from vplants.self_similarity.cst import *
from vplants.self_similarity.mtg2graph import *"""

from vplants.tree_reduction.graph import *
from openalea.container.graph import *
import string

"""mtg_filename = './example_MTG.mtg'

g = MTG(mtg_filename)

print 'Create the graph'
tree = mtg2graph(1, 2)

print 'Reduction'
dag_exact = tree_reduction(tree)
tree_rec = tree_reconstruction(tree)"""

"""   rec directe
def test1():
    edges=[(1,2),(2,3),(3,4),(3,5),(4,6),(4,5),(5,6)]
    edges_weight={0:1,1:2,2:1,3:1,4:1,5:1,6:1}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 6
    assert g.nb_edges() == 7
    return g"""

""" rec indirecte"""
def test1():
    edges=[(1,2),(1,7),(2,3),(2,6),(3,4),(3,7),(4,6),(4,5),(5,7),(6,8),(7,8)]
    edges_weight={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:3,10:2}
    root = 1
    g = from_edges(root, edges, edges_weight)
    assert g.nb_vertices() == 8
    assert g.nb_edges() == 11
    return g

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
    QIP=[]
    for eid in Quasi_isom_desc_vertices(vid, g, signat_edges):
        QIP.append([eid])
    while (len(QIP)>1):
        if len(QIP)==0:
            QPP=QIP[:]
        else:
            QPP=[]
            QPP.append(QIP[0][:])
            Q=QIP[:]
            a=True
            while a==True:
                elt=QPP[len(QPP)-1]
                eid=elt[len(elt)-1]
                i=0
                while i<len(Q):
                    eid2=Q[i][0]
                    if g._edges[eid][1]==g._edges[eid2][0]:
                        QPP.append(Q[i])
                        Q=Q[0:i]+Q[i+1:]
                        a=True
                    else:
                        i=i+1
                        a=False
            if len(QPP)>1:
                M.append(QPP)

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

def MQPPRT(g, signat_edges):
    '''
    '''
    MRT=[]
    for vid in g.vertices():
        MRT=MRT+MQPP(vid, g, signat_edges)
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

def equal_string(a,b):
    if (a=='alpha' and b=='alpha'):
        e=False
    elif a=='' and b=='':
        e=True
    elif (a=='' and b!='') or (a!='' and b==''):
        e=False
    else:
        sub=a[0:4]
        nsub=string.count(a,sub)
        if nsub!=string.count(b,sub):
            e=False
        else:
            a=string.replace(a,sub,'')
            b=string.replace(b,sub,'')
            e=equal_string(a,b)
    return e

def sup_edge(x,y,list_edge,list_node,G_edge,nbr_edge):
    """supprime l'arc entre les noeuds x et y"""
    k=ind_find(x,list_edge)
    if y in list_edge[k][1]:
        list_edge[k][1].remove(y)
        if (nbr_edge[G_edge[(list_node[x],list_node[y])]]>=1):
            nbr_edge[G_edge[(list_node[x],list_node[y])]]-=1
        else:
            nbr_edge.pop(G_edge[(list_node[x],list_node[y])])
        G_edge.pop((list_node[x],list_node[y]))
    return list_edge,G_edge,nbr_edge

def nbr_parent(i,list_edge):
    """retourne le nombre des noeuds parents du noeud i"""
    k=0
    nbr=0
    while k<len(list_edge):
        if i in list_edge[k][1]:
            nbr=nbr+1
        k=k+1
    return nbr

def list_parent(i,list_edge):
    """retourne le nombre des noeuds parents du noeud i"""
    k=0
    lp=[]
    while k<len(list_edge):
        if i in list_edge[k][1]:
            lp.append(list_edge[k][0])
        k=k+1
    return lp
        
def ind_find(x,list_edge):
    """recherche l'indice dans list_edge de l'element associé au noeud x
    (pour accéder a la liste de ces fils)"""
    k=0
    while (list_edge[k][0]!=x) and (k<len(list_edge)):
            k=k+1  
    return k

def ajust_sup(list_edge,G_edge,list_node,G_node,nbr_edge):
    """Après la suppression d'un arc le graphe G a besoin d'un ajustement
    en supprimant tous les nœuds inaccessibles qui ne possède aucun parent"""
    k=1
    mod=True
    while mod==True:
        mod=False
        while (k<len(list_edge)):
            if (nbr_parent(list_edge[k][0],list_edge)==0):
                mod=True
                G_node.pop(list_node[list_edge[k][0]])
                k1=0
                while (k1<len(list_edge[k][1])):
                    if (nbr_edge[G_edge[(list_node[list_edge[k][0]],list_node[list_edge[k][1][k1]])]]>=1):
                        nbr_edge[G_edge[(list_node[list_edge[k][0]],list_node[list_edge[k][1][k1]])]]-=1
                    else:
                        nbr_edge.pop(G_edge[(list_node[list_edge[k][0]],list_node[list_edge[k][1][k1]])])
                    G_edge.pop((list_node[list_edge[k][0]],list_node[list_edge[k][1][k1]]))
                    k1+=1
                list_node=list_node[0:list_edge[k][0]]+list_node[list_edge[k][0]+1:len(list_node)]
                for k1 in range(len(list_edge)):
                    if list_edge[k1][0]>list_edge[k][0]:
                        list_edge[k1][0]-=1
                    for k2 in range(len(list_edge[k1][1])):
                        if list_edge[k1][1][k2]>list_edge[k][0]:
                            list_edge[k1][1][k2]-=1
                list_edge=list_edge[0:k]+list_edge[k+1:len(list_edge)]
            else: k+=1
    return list_edge,G_edge,list_node,G_node,nbr_edge


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
  
def creer_symbol(V): #ajouter un nouveau symbole correspondant a un noeud
  c=''
  for c in string.lowercase+string.uppercase:
    if c not in V and c!='[' and c!='F' and c!=']' and c!='l' and c!='W':
      break
  return c
        
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
    
def nbr_arc_construction(G_edge):
    nbr_edge={}
    d=G_edge.values()
    while d!=[]:
        val=d[0]
        d.remove(val)
        nbr=1
        k=0
        while (k<len(d)):
          if equal_string(val,d[k]):
            nbr=nbr+1
            d.remove(d[k])
          else: k=k+1  
        nbr_edge[val]=nbr
    nbr_edge['alpha']=1
    return nbr_edge


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
                print '********'
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

   
