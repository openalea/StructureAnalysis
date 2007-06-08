/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (aida.ouangraoua@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/Ordered/matching_GOQDE.cpp,v $
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */


#include "matching_GOQDE.h"

//const int EMPTY_TREE=-1;
//const int EMPTY_NODE=-1;
const int EMPTY_CLASS=-1;


//-----------------------------------------------------------------------------
// Classe permettant le calcul d'une valeur minimale parmi plusieurs 
//-----------------------------------------------------------------------------


class Min {
  DistanceType value;
public:
  Min(DistanceType v):value(v) {}
  Min &operator<<(DistanceType v) {
    if (v<value) value=v;
    return *this;
  }
  operator DistanceType() {
    return value;
  }
};



// -------------
  // Constructeur
  // -------------
MatchingGoqde::MatchingGoqde(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  Matching_O(input,reference,nodeDistance)
{
  T1=&input;
  T2=&reference;
  CostMatrix=&nodeDistance;
  int sizeT1 = T1->getNbVertex();
  int sizeT2 = T2->getNbVertex();

 //tableau de distance entre forêts
   int sizeFD = sizeT1*sizeT2;   
   fDist.resize(sizeFD);
   fDistC1.resize(sizeFD);
   fDistC2.resize(sizeFD);
   fDistC1C2.resize(sizeFD);
   
   L1.resize(sizeT1);
   L2.resize(sizeT2);
  //sequence correspondant à l'alignement
   alignement = new Sequence(); 

//Coûts associés à l'alignement trouvé
   insertCost = 0;
   deleteCost = 0;
   matchCost = 0; 

}
// -------------

// Destructeur
  // -------------
MatchingGoqde::~MatchingGoqde()
{
  
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les forets F1[L1[input_vertex]..i1] et F2[L2[reference_vertex]..j1]
// ----------------------------------------------------------------------------------


void MatchingGoqde::distanceBetweenForest(int input_vertex,int reference_vertex)
{
    int i,j,i1,j1;
   
    i= input_vertex;
    j= reference_vertex;   
    //Calcul récursif de la distance entre les forêts T1[L1[i]..i] et T2[L2[j]..j]  
    int fNum = fix(i,j);
    fDist[fNum][0][0] = 0;
    
    //cout<<"....Computing Delete Cost"<<endl;
    for (i1 = L1[i]; i1 >= i; i1--) {
      fDist[fNum][idx(L1[i],i1)][0] = fDist[fNum][idx(L1[i],i1+1)][0] + getDeleteCost(i1);
    }
    
    //cout<<"....Computing Insert Cost"<<endl;
    for (j1 = L2[j]; j1 >= j; j1--) {
      fDist[fNum][0][idx(L2[j],j1)] = fDist[fNum][0][idx(L2[j],j1+1)] + getInsertCost(j1);
    }
    
    for (i1 = L1[i]; i1 >= i; i1--) {
      for (j1 = L2[j]; j1 >= j; j1--) {
	//cout<<"....Computing : "<<i1<<" "<<j1<<endl;
	//cout<<"....-1 -1"<<endl;
	computeDist(L1[i],i1, L2[j],j1,-1,-1);
	//cout<<"....C1 -1"<<endl;	
	computeDist(L1[i],i1, L2[j],j1,classRoot(1,lca(1,L1[i],i1)),-1);
	//cout<<"....-1 C2"<<endl;	
	//cout<<"....lca : "<<lca(2,L2[j],j1)<<endl;
	//cout<<"....root : "<<classRoot(2,lca(2,L2[j],j1))<<endl;

	computeDist(L1[i],i1, L2[j],j1,-1,classRoot(2,lca(2,L2[j],j1)));
	//cout<<"....C1 C2"<<endl;	
	computeDist(L1[i],i1, L2[j],j1,classRoot(1,lca(1,L1[i],i1)),classRoot(2,lca(2,L2[j],j1)));	
      }
    }
}

void MatchingGoqde::computeDist(int li, int i, int lj, int j, int c1, int c2)
{
  
  DistanceType match,dele,inse,min;
  int m1,m2;
  int fNum = fix(getKey1(li),getKey2(lj));
  
  m1 = Min(li)<<L1[classRoot(1,i)]; 
  //cout<<"........m1 : "<<m1;
  
  m2 = Min(lj)<<L2[classRoot(2,j)];  
  //cout<<" m2 : "<<m2<<endl; 
  
  //Calcul de FDist sans contrainte
  match =  getDist(li,m1+1,lj,m2+1, c1, c2) + getDist(m1,L1[i]+1,m2,L2[j]+1,classRoot(1,i),classRoot(2,j)) + getDist(L1[i],i+1,L2[j],j+1,classRoot(1,i),classRoot(2,j)) +  getSubstitutionCost(i,j);
  //cout<<"........match : "<<match<<endl; 
  

  dele = getDist(li,i+1,lj,j,c1,c2) +  getDeleteCost(i);
  //cout<<"........dele : "<<dele<<endl; 
  
  inse = getDist(li,i,lj,j+1,c1,c2) +  getInsertCost(j);
  //cout<<"........inse : "<<inse<<endl; 
  
  if((classRoot(1,i) == c1 && classRoot(2,j) == c2) || (classRoot(1,i) != c1 && classRoot(2,j) != c2))
    min = Min(dele)<<inse<<match;
  else
    min = Min(dele)<<inse;
  
  if(c1 == -1 && c2 == -1)
    fDist[fNum][idx(li,i)][idx(lj,j)] = min;
  if(c1 != -1 && c2 == -1)
    fDistC1[fNum][idx(li,i)][idx(lj,j)] = min;
  if(c1 == -1 && c2 != -1)
    fDistC2[fNum][idx(li,i)][idx(lj,j)] = min;
  if(c1 != -1 && c2 != -1)
    fDistC1C2[fNum][idx(li,i)][idx(lj,j)] = min;
  //cout<<"....End ComputeDist"<<endl;
}


//Operator
//Distance contrainte entre deux arbres
DistanceType MatchingGoqde::getDBT(int input_vertex,int reference_vertex,int c1, int c2) {
  return   getDist(L1[input_vertex],input_vertex, L2[reference_vertex],reference_vertex,c1,c2);
} 



//Distance contrainte entre deux forets
DistanceType MatchingGoqde::getDist(int input_vertex1,int input_vertex2,int reference_vertex1,int reference_vertex2, int c1, int c2) {
  return getFromFDist(input_vertex1,input_vertex2,reference_vertex1,reference_vertex2, c1,c2);
  
} 


DistanceType MatchingGoqde::getFromFDist(int input_vertex1,int input_vertex2,int reference_vertex1,int reference_vertex2, int c1, int c2) {
  DistanceType dist;
  
  if ((input_vertex1 >= input_vertex2) && (reference_vertex1 >= reference_vertex2)){
    
    if(c1 != classRoot(1,lca(1,input_vertex1,input_vertex2)))
      c1 = -1;
    
    if(c2 != classRoot(2,lca(2,reference_vertex1,reference_vertex2)))
      c2 = -1;
    
    if(c1==-1 && c2 == -1)
      dist = fDist[fix(getKey1(input_vertex1),getKey2(reference_vertex1))][idx(input_vertex1,input_vertex2)][idx(reference_vertex1,reference_vertex2)]; 
    if(c1!=-1 && c2 == -1)
      dist = fDistC1[fix(getKey1(input_vertex1),getKey2(reference_vertex1))][idx(input_vertex1,input_vertex2)][idx(reference_vertex1,reference_vertex2)]; 
    if(c1==-1 && c2 != -1)
      dist = fDistC2[fix(getKey1(input_vertex1),getKey2(reference_vertex1))][idx(input_vertex1,input_vertex2)][idx(reference_vertex1,reference_vertex2)]; 
    if(c1!=-1 && c2 != -1)
      dist = fDistC1C2[fix(getKey1(input_vertex1),getKey2(reference_vertex1))][idx(input_vertex1,input_vertex2)][idx(reference_vertex1,reference_vertex2)]; 
  }
  
  if((input_vertex1 < input_vertex2)  && (reference_vertex1 >= reference_vertex2))
    dist =  fDist[fix(getKey1(input_vertex1), getKey2(reference_vertex1))][0][idx(reference_vertex1,reference_vertex2)];
  
  if ((input_vertex1 >= input_vertex2) && (reference_vertex1 < reference_vertex2))
    dist = fDist[fix(getKey1(input_vertex1),getKey2(reference_vertex1))][idx(input_vertex1,input_vertex2)][0]; 
  
  if ((input_vertex1 < input_vertex2) && (reference_vertex1 < reference_vertex2))
    dist = 0;
  
  return dist;  
}


//Calcul des feuilles gauches associées à chaque sommet dans un arbre
void MatchingGoqde::computeLeft(int numT, int root){
  int leaf;
  int node;
  TreeGraph* T;
  if(numT == 1)
    T = T1;
  else
    T = T2;
  int nbchild = T->getNbChild(root);
  
  if (nbchild == 0)
    {
      //si root est une feuille, on s'arrête là
      leaf = root;
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
      
      int leftchild = T->child(root,nbchild);
      if(numT == 1)
	leaf = L1[leftchild];
      else
	leaf = L2[leftchild];
    }   
  if(numT == 1)
    L1[root]=leaf;
  else
    L2[root]=leaf;
} 



//Calcul des numéros des racines particulières dans un arbre
void MatchingGoqde::computeKeyroots(int numT){
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
  
  for(int i=Lsize-1;i>=0;i--){
    //Pour chaque noeud  de T numéroté i suivant la numérotation postfixe
    father = T->father(i);
    if(father != -1){
      //si i n'est pas la racine
      //si L[père[i]] != L[i] alors i est un keyroot 
      if(numT == 1){   
	if(L1[father] != L1[i]){
	  keyroots1.resize(j+1);
	  keyroots1[j]=i;
	  j++;    
	}    
      }
      else{
	if(L2[father] != L2[i]){
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




//Calculs init (leftleaves et keyroots)
void MatchingGoqde::precompute(){
  
  int i,j;
  int sizeClass1,sizeClass2;
  
  computeLeft(1,0);
  
  /*cout<<"L1 : ";
  for (i = 0; i < T1->getNbVertex(); i++){
    cout<<" "<<L1[i];
  }
  cout<<endl;
  */
  computeKeyroots(1);
  sizeKeyroots1= keyroots1.size();
  /*
  cout<<"keyroots1 : ";
  for (i = 0; i < sizeKeyroots1; i++){
    cout<<" "<<keyroots1[i];
  }
  cout<<endl;
  */

  computeLeft(2,0);
  /*
  cout<<"L2 : ";
  for (i = 0; i < T2->getNbVertex(); i++){
    cout<<" "<<L2[i];
  }
  cout<<endl;
  */
  computeKeyroots(2);
  sizeKeyroots2= keyroots2.size();  
  /*
  cout<<"keyroots2 : ";
  for (i = 0; i < sizeKeyroots2; i++){
    cout<<" "<<keyroots2[i];
  }
  cout<<endl;
  */
  
  for(  i = 0; i < sizeKeyroots1; i++){
    for( j = 0; j < sizeKeyroots2; j++){ 
      fDist[fix(keyroots1[i],keyroots2[j])].resize(L1[keyroots1[i]]-keyroots1[i]+2);
      fDistC1[fix(keyroots1[i],keyroots2[j])].resize(L1[keyroots1[i]]-keyroots1[i]+2);
      fDistC2[fix(keyroots1[i],keyroots2[j])].resize(L1[keyroots1[i]]-keyroots1[i]+2);
      fDistC1C2[fix(keyroots1[i],keyroots2[j])].resize(L1[keyroots1[i]]-keyroots1[i]+2);
      
      for(int k = 0; k < L1[keyroots1[i]]-keyroots1[i]+2; k++){ 
	fDist[fix(keyroots1[i],keyroots2[j])][k].resize(L2[keyroots2[j]]-keyroots2[j]+2);
	fDistC1[fix(keyroots1[i],keyroots2[j])][k].resize(L2[keyroots2[j]]-keyroots2[j]+2);
	fDistC2[fix(keyroots1[i],keyroots2[j])][k].resize(L2[keyroots2[j]]-keyroots2[j]+2);
	fDistC1C2[fix(keyroots1[i],keyroots2[j])][k].resize(L2[keyroots2[j]]-keyroots2[j]+2);
      }
    }
  }   
  //cout<<"Precompute End"<<endl;  
} 


int MatchingGoqde::fix(int a, int b){
 int sizeT2 = T2->getNbVertex(); 
  return a*sizeT2 +b;
}



int MatchingGoqde::idx(int a , int b){
return (a<b)?0:(a-b+1);
}

int MatchingGoqde::getKey1(int li){
  int i;
  for (i = sizeKeyroots1 - 1; L1[keyroots1[i]] != L1[li]; i--);
  return keyroots1[i];
}

int MatchingGoqde::getKey2(int lj){
  int j;
  for (j = sizeKeyroots2 - 1; L2[keyroots2[j]] != L2[lj]; j--);
  return keyroots2[j];
}

int  MatchingGoqde::lca(int numT, int i1, int i2){
  
  int father;
  int n;
  
  n = i2;
  if(numT == 1){
    while (L1[n]!=L1[i1])
      {
	n = T1->father(n);
      }
  }
  else{
    while (L2[n]!=L2[i1])
      {
	n = T2->father(n);
      }
  }
  return n;
}

int  MatchingGoqde::classRoot(int numT, int i){
  
  int father;
  int r;
  
  
  father = i;
  if(numT == 1){
    while (T1->getNode(father)->getComplex()==T1->getNode(i)->getComplex() && father != 0)
      {
	r = father;
	father = T1->father(father);
      }
    if(T1->getNode(father)->getComplex()==T1->getNode(i)->getComplex())
      r = father;
  }
  else{
    
    while (T2->getNode(father)->getComplex()==T2->getNode(i)->getComplex() && father != 0)
      {
	r = father;
	father = T2->father(father);
      }
    if(T2->getNode(father)->getComplex()==T2->getNode(i)->getComplex())
      r = father;
  }
  
  return r;
}


DistanceType  MatchingGoqde::match(){

  DistanceType dist;
  
  int i,j,v,w;
  int size1,size2;
  
  precompute();


  size1 = sizeKeyroots1;
  size2 = sizeKeyroots2;
  //cout<<"Avant Boucle"<<endl;
  //Algo de Zhang et Shasha
  for (i = 0; i < sizeKeyroots1; i++){
    for (j = 0; j < sizeKeyroots2; j++){
      v =  keyroots1[i];
      w =  keyroots2[j];
      //cout<<"Computing : "<<v<<" "<<w<<endl;
      distanceBetweenForest(v, w);    
    }
  }
  //cout<<"Après Boucle"<<endl;
  v =  keyroots1[sizeKeyroots1-1];
  w =  keyroots2[sizeKeyroots2-1];
  dist = getDBT(v,w,-1,-1);
 
  cout <<"Distance : "<<dist<<endl;
  return dist;
} 

