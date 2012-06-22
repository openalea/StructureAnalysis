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
