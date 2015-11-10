#include "MS_O_Matching.h"
#include <iostream>
#include <strstream>
#include <cstdio>
#include <cmath>


// ----------------------------------------------------------------------------------
// Calcul des feuilles les plus à gauches (l(k)) pour tous les noeuds
// de l'arborescence de racine root  
// ------------------------------------------------------------------------


void MS_O_Matching::computeLeft(int numT, int root)
{
  int leaf;
  int node;
  TreeNode* nodeSon;
  TreeGraph* T;
  if(numT == 1)
    T = T1;
  else
    T = T2;
  int nbchild = T->getNbChild(root);
  TreeNode* _node=T->getNode(root);

  if (nbchild == 0)
    {
      //si root est une feuille, on s'arrête là
      leaf = _node->getNumPostfix();
    }
  else
    {
      //sinon lancement de computeLeft pour tous les fils de root
      // et L[root] = L[fils_gauche_de_root]           
      int child =1;
      while (child <= nbchild) {
	node = T->child(root,child);
	computeLeft(numT,node);
	child++;
      }
      
      int leftchild = T->child(root,1);
      nodeSon = T->getNode(leftchild);
      if(numT == 1)
	leaf = L1[nodeSon->getNumPostfix()];
      else
	leaf = L2[nodeSon->getNumPostfix()];
    }   
  if(numT == 1)
    L1[_node->getNumPostfix()]=leaf;
  else
    L2[_node->getNumPostfix()]=leaf;
}


// ----------------------------------------------------------------------------------
// Calcul des racines particulières (keyroots) 
// ----------------------------------------------------------------------------------


void MS_O_Matching::computeKeyroots(int numT)
{

  int father;
  int Lsize;
  TreeGraph* T;

  if(numT == 1){
    Lsize = L1.size();
    T = T1;
  }
  else{
    Lsize = L2.size();
    T = T2;
  }
  
  int j =0;
  
  for(int i=0;i<=Lsize-1;i++){
    //Pour chaque noeud  de T numéroté i suivant la numérotation postfixe
    int num = T->getRealNumber(i);
    father = T->father(num);
    if(father != -1){
      //si i n'est pas la racine
      //si L[père[i]] != [i] alors i est un keyroot 
      TreeNode* node = T->getNode(father);
      if(numT == 1){   
	if(L1[node->getNumPostfix()] != L1[i]){
	  keyroots1.resize(j+1);
	  keyroots1[j]=i;
	  j++;    
	}    
      }
      else{
	if(L2[node->getNumPostfix()] != L2[i]){
	  keyroots2.resize(j+1);
	  keyroots2[j]=i;
	  j++;    
	}    
      }
    }
    else       
      {
	//si i est la racine alors i est keyroot
	if(numT == 1){   
	  keyroots1.resize(j+1);
	  keyroots1[j]=i;
	}
	else{
	  keyroots2.resize(j+1);
	  keyroots2[j]=i;
	}
	j++;  
      }
  }
}

//--------------------------------------------------------------------------------
// Précalcul des l(k) et keyroots avant lancement de l'algorithme
//------------------------------------------------------------------------
void MS_O_Matching::precompute()
{
  int i,j;
  int sizeT1, sizeT2;
  int root1, root2;
  
  
  sizeT1 = T1->getNbVertex();
  sizeT2 = T2->getNbVertex();
  
  
  L1.resize(sizeT1);
  computeLeft(1,0);
  
  computeKeyroots(1);
  sizeKeyroots1= keyroots1.size();
  
  
  L2.resize(sizeT2);
  computeLeft(2,0);
  
  computeKeyroots(2);
  sizeKeyroots2= keyroots2.size();
 
  for(  i = 0; i < sizeKeyroots1; i++){
    for( j = 0; j < sizeKeyroots2; j++){ 
      //cout<<keyroots1[i]<<";"<<keyroots2[j]<<endl;
      fDist[fix(keyroots1[i],keyroots2[j])].resize(keyroots1[i]-L1[keyroots1[i]]+2);
      
      for(int k = 0; k < keyroots1[i]-L1[keyroots1[i]]+2; k++){ 
	fDist[fix(keyroots1[i],keyroots2[j])][k].resize(keyroots2[j]-L2[keyroots2[j]]+2);
      }
    }
  }
}
