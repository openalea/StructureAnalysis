# -*- coding: cp1252 -*-
import string

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


def DAG_init():   
    node=[1,2,3,4,5,6,7,8,9,10,11]
    edge={1:[2,8],2:[3,4],3:[11],4:[5,6],5:[11],6:[3,7],7:[5,8],8:[3,9],9:[5,10],10:[11]}
    val={(1,2):1,(1,8):1,(2,4):1,(2,3):1,(3,11):2,(4,5):1,(4,6):1,(5,11):3,(6,7):1,(6,3):1,(7,5):1,(7,8):1,(8,3):1,(8,9):1,(9,5):1,(9,10):1,(10,11):1}
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
        if functions.equal_string(G_edge[G_edge.keys()[k]],nbr_edge.keys()[j]):
          G_edge[G_edge.keys()[k]]=nbr_edge.keys()[j]    
    return list_node,list_edge,G_node,G_edge,nbr_edge,parameters,parameters1
    
def nbr_arc_construction(G_edge):
    nbr_edge={}
    d=G_edge.values()
    while d!=[]:
        val=d[0]
        d.remove(val)
        nbr=1
        k=0
        while (k<len(d)):
          if functions.equal_string(val,d[k]):
            nbr=nbr+1
            d.remove(d[k])
          else: k=k+1  
        nbr_edge[val]=nbr
    nbr_edge['alpha']=1
    return nbr_edge

