#include "MS_O_Matching.h"
#include <iostream>
#include <cstdio>
#include <cmath>

//--------------------------------------------------------------------------------
// Calcul du coût local d'insertion d'un noeud
//------------------------------------------------------------------------
DistanceType MS_O_Matching::insLocalNode(TreeGraph* T, int w){
  DistanceType Dist;
  
  int wReal = T->getRealNumber(w);  
  TreeNode* _node =  T->getNode(wReal);
  
  //Si il ne s'agit pas de la dernière échelle de comparaison 
  if(T->getMTG()->vscale(_node->getVertex()) < _init){
    VId root = _node->getVertex();
    TreeGraph* Tree2 = new TreeGraph(*(T->getMTG()),root,COMPO);
    int nbVertex = Tree2->getNbVertex();
    
    //somme des coûts d'insertion des sous-noeuds à l'échelle courante+1 
    Dist = 0 ;
    for(int i =0; i<nbVertex;i++ ){    
      Dist = Dist + insLocalNode(Tree2,i);
    } 
    
    delete (TreeGraph*) Tree2;
    
  }
  else
    {
      //sinon coût d'insertion suivant la fonction coût utilisé
      Dist = CostMatrix->getInsertionCost(_node);
    }
  
  
  return Dist;
}

//--------------------------------------------------------------------------------
// Calcul du coût local de suppression d'un noeud
//------------------------------------------------------------------------
DistanceType MS_O_Matching::delLocalNode(TreeGraph* T, int v){
  DistanceType Dist;
  
  
  int vReal = T->getRealNumber(v);
  TreeNode* _node =  T->getNode(vReal);
  
  //Si il ne s'agit pas de la dernière échelle de comparaison 
  if(T->getMTG()->vscale(_node->getVertex()) < _init){
    VId root = _node->getVertex();
    TreeGraph* Tree1 = new TreeGraph(*(T->getMTG()),root,COMPO);
    int nbVertex = Tree1->getNbVertex();
    
    //somme des coûts de suppression des sous-noeuds à l'échelle courante+1 
    Dist = 0 ;
    for(int i =0; i<nbVertex;i++ ){
      Dist = Dist + delLocalNode(Tree1,i);
    } 
    
    delete (TreeGraph*) Tree1;
    
  } 
  else
    {      
      //sinon coût de suppression suivant la fonction coût utilisé
      Dist = CostMatrix->getDeletionCost(_node);
    }
  
  return Dist;
}

//--------------------------------------------------------------------------------
// Calcul du coût local de remplacement d'un noeud par un autre
//------------------------------------------------------------------------
DistanceType MS_O_Matching::matchLocalNode(int v, int w){
  DistanceType Dist;
  int vReal = T1->getRealNumber(v);
  int wReal = T2->getRealNumber(w);
  
  TreeNode* _node1 =  T1->getNode(vReal);
  TreeNode* _node2 =  T2->getNode(wReal);
  
  //Si il ne s'agit pas de la dernière échelle de comparaison 
  if(T1->getMTG()->vscale(_node1->getVertex()) < _init){
    VId root1 = _node1->getVertex();
    TreeGraph* Tree1 = new TreeGraph(*(T1->getMTG()),root1,COMPO);
    
    
    VId root2 = _node2->getVertex();
    TreeGraph* Tree2 = new TreeGraph(*(T2->getMTG()),root2,COMPO);
    
    
    //distance entre les sous-arborescences à l'échelle courante+1 
    MS_O_Matching* M=new MS_O_Matching(*Tree1,*Tree2,*CostMatrix,_init);
    Dist=M->match();
    delete (MS_O_Matching*) M;
    delete (TreeGraph*) Tree1;
    delete (TreeGraph*) Tree2;
  }
  else
    {
      //sinon coût de remplacement suivant la fonction coût utilisé
      Dist =  CostMatrix->getChangingCost(_node1,_node2);
    }
  return Dist;
}
