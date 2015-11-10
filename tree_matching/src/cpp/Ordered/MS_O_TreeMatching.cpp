#include "MS_O_Matching.h"
#include <iostream>
#include <cstdio>
#include <cmath>




// -------------
// Constructeur
// -------------
MS_O_Matching::MS_O_Matching(TreeGraph& input,TreeGraph& reference, NodeCost& nodeDistance, int init)
{
  T1=&input;
  T2=&reference;
  
  int sizeT1 = T1->getNbVertex();
  int sizeT2 = T2->getNbVertex();
  
  //Fonction de coût  
  CostMatrix = &nodeDistance;

  //flag indiquant le premier objet MS_O_Matching
  _init = init;

  //tableau de distance entre forêts
   int sizeFD = sizeT1*sizeT2;   
   fDist.resize(sizeFD);
   
   //sequence correspondant à l'alignement
   alignement = new Sequence(); 

   //Coûts associés à l'alignement trouvé
   insertCost = 0;
   deleteCost = 0;
   matchCost = 0; 
   
   symbols2["N"] = 0;
   symbols2["I"] = 1;
   symbols2["B"] = 2;
   symbols2["H"] = 3;
   symbols2["M"] = 4;
   symbols2["S"] = 5;
   symbols2["-"] = 6;

   symbols["-"] = 0;
   symbols["A"] = 1;
   symbols["U"] = 2;
   symbols["G"] = 3;
   symbols["C"] = 4;
   symbols["AU"] = 5;
   symbols["UA"] = 6;
   symbols["AG"] = 7;
   symbols["GA"] = 8;
   symbols["AC"] = 9;
   symbols["CA"] = 10;
   symbols["UG"] = 11;
   symbols["GU"] = 12;
   symbols["UC"] = 13;
   symbols["CU"] = 14;
   symbols["GC"] = 15;
   symbols["CG"] = 16;
   symbols["AA"] = 17;
   symbols["UU"] = 18;
   symbols["GG"] = 19;
   symbols["CC"] = 20;
   symbols["N"] = 21;
}

// -------------
// Destructeur
// -------------
MS_O_Matching::~MS_O_Matching()
{
}

const DistanceType MS_O_Matching::epsilon = 0.01;

DistanceType MS_O_Matching::CostMatrix2[8][8]
 = {
   {0,10000,10000,10000,10000,10000,10000,10000},
   {10000,0,3,8,8,10000,5},
   {10000,3,0,8,8,10000,5},
   {10000,8,8,0,8,10000,100},
   {10000,8,8,8,0,10000,75},
   {10000,10000,10000,10000,10000,0,5},
   {10000,5,5,100,75,5,0},
 };


DistanceType MS_O_Matching::CostMatrix1[NB_SYMBOLS][NB_SYMBOLS]
 = {
  {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
  {1, 0, 1, 1, 1, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
  {1, 1, 0, 1, 1, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
  {1, 1, 1, 0, 1, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
  {1, 1, 1, 1, 0, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000},
  {1, 10000, 10000, 10000, 10000, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 10000},
  {1, 10000, 10000, 10000, 10000, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 10000},
  {0, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 0},

};


/*int MS_O_Matching::CostMatrix1[NB_SYMBOLS][NB_SYMBOLS]
 = {
  {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1},
  {1, 0, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {1, 2, 0, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {1, 2, 2, 0, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {1, 2, 2, 2, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3},
  {1, 3, 3, 3, 3, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2},
  {1, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0},

};
*/




// ----------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
void MS_O_Matching::computeTreeDist(int input_vertex,int reference_vertex)
{

  //Calcul récursive des distances entre les sous-forêts
  computeForestsDistances(input_vertex, reference_vertex);
  //int d = fDist[fix(getKey1(input_vertex), getKey2(reference_vertex))][idx(L1[input_vertex],input_vertex)][idx(L2[reference_vertex],reference_vertex)]; 
}

//------------------------------------------------------------------------------------
//Distance entre deux arbres indexés par leurs racines
//-------------------------------------------------------------------------------------
DistanceType MS_O_Matching::getTreeDistance(int input_vertex,int reference_vertex)
{
  return   fDist[fix(getKey1(input_vertex), getKey2(reference_vertex))][idx(L1[input_vertex],input_vertex)][idx(L2[reference_vertex],reference_vertex)]; 
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  MS_O_Matching::match()
{
  DistanceType dist;
  
  int i,j,v,w;
  int size1,size2;
  
  
  precompute();
  
  size1 =sizeKeyroots1;
  size2 =sizeKeyroots2;
  
  //Algo de Zhang et Shasha
  for (i = 0; i < sizeKeyroots1; i++){
    for (j = 0; j < sizeKeyroots2; j++){
      v =  keyroots1[i];
      w =  keyroots2[j];
      computeTreeDist(v, w);
    }
  }
  
  v =  keyroots1[sizeKeyroots1-1];
  w =  keyroots2[sizeKeyroots2-1];
  dist = getTreeDistance(v,w);
  
  //Calcul d'un alignement optimal
  makeAlignement(keyroots1[size1-1],keyroots2[size2-1]);
  
  
  return dist;
}


//-----------------------------------------------------------
// Calcul d'un alignement optimal entre les arborescences de 
// racines étiquetés i et j suivant la numérotation postfixe
//------------------------------------------------------------ 

void MS_O_Matching::makeAlignement(int i, int j){
  
  int li = L1[i];
  int lj = L2[j];
  subAlign(li, i, lj, j);
}

//-----------------------------------------------------------
// Calcul d'un alignement optimal entre les forëts de 
// T1[li..i] et T2[lj..j] suivant la numérotation postfixe
//------------------------------------------------------------ 

void MS_O_Matching::subAlign(int li, int i, int lj, int j){
 
  if (L1[i] < li || L2[j] < lj) {
    return;
  }
  
  if(i==-1){
    for(int k=lj; k<=j; k++){
      insertCost +=insLocalNode(T2,k); 
      alignement->putNbIns(alignement->getNbIns()+getNbVertexAtScaleInit(T2, k)); 
    }
    return;
  }
  
  if(j==-1){
    for(int k=li; k<=i; k++){
      deleteCost +=delLocalNode(T1,k);
      alignement->putNbDel(alignement->getNbDel()+getNbVertexAtScaleInit(T1, k));  
    }
    return;
  }
  
  //Récupération des numéros préfixes associés aux noeuds
  int iReal = T1->getRealNumber(i);
  int jReal = T2->getRealNumber(j);
  TreeNode* _node1 =  T1->getNode(iReal);
  TreeNode* _node2 =  T2->getNode(jReal);
  
  if(li == i && lj == j){  
    if(T1->getMTG()->vscale(_node1->getVertex()) == _init){
      //Si il s'agit de la dernière échelle de comparaison, on ajoute le couple de noeuds      
      alignement->append( _node1->getNumber(),_node2->getNumber(), matchLocalNode(i,j));   
      matchCost += matchLocalNode(i,j);
    }
    else
      {
	//sinon il faut ajouter l'alignement correspondant aux sous-arborescences
	//des noeuds i et j à l'échelle courante+1.
	
	VId root1 = _node1->getVertex();
	TreeGraph* Tree1 = new TreeGraph(*(T1->getMTG()),root1,COMPO);
	
	
	VId root2 = _node2->getVertex();
	TreeGraph* Tree2 = new TreeGraph(*(T2->getMTG()),root2,COMPO);
	
	MS_O_Matching* M=new MS_O_Matching(*Tree1,*Tree2,*CostMatrix,_init);
	M->match();
	alignement->add(M->getSequence());
	matchCost += M->getMatchingCost();
	insertCost += M->getInsertionCost();
	deleteCost += M->getDeletionCost();
	alignement->putNbIns(alignement->getNbIns()+M->getSequence()->getNbIns());
	alignement->putNbDel(alignement->getNbDel()+M->getSequence()->getNbDel());
	delete (MS_O_Matching*) M;
	delete (TreeGraph*) Tree1;
	delete (TreeGraph*) Tree2;
      }
    return;
  } 
  
  
  int fNum = fix(getKey1(li),getKey2(lj));
  DistanceType f =  fDist[fNum][idx(li,i)][idx(lj,j)];
  DistanceType fi_1j_1 = fDist[fNum][idx(li,i-1)][idx(lj,j-1)];
  DistanceType fi_1 = fDist[fNum][idx(li,i-1)][idx(lj,j)];
  DistanceType fj_1 = fDist[fNum][idx(li,i)][idx(lj,j-1)];
  
  
  if(L1[i]==li && L2[j]==lj) {
    // si foret[li...i] et foret[lj...j] sont des arbres 
    
    if (f - fi_1 - delLocalNode(T1,i) == 0) {
      subAlign(li,i-1,lj,j);
      deleteCost+=delLocalNode(T1,i);      
      alignement->putNbDel(alignement->getNbDel()+getNbVertexAtScaleInit(T1, i)); 
      return ;
    } 
    
    if (f - fj_1 - insLocalNode(T2,j) == 0) {
      subAlign(li,i,lj,j-1);
      insertCost +=insLocalNode(T2,j);
      alignement->putNbIns(alignement->getNbIns()+getNbVertexAtScaleInit(T2, j)); 
      return;
    } 
    
    if (f - fi_1j_1 - matchLocalNode(i,j) == 0) {
      
      if(T1->getMTG()->vscale(_node1->getVertex()) == _init){
	//Si il s'agit de la dernière échelle de comparaison, on ajoute le couple de noeuds
	//cout<<matchLocalNode(i,j)<<endl;
	alignement->append( _node1->getNumber(),_node2->getNumber(), matchLocalNode(i,j));
	matchCost+= matchLocalNode(i,j);  
      }
      else
	{
	  //sinon il faut ajouter l'alignement correspondant aux sous-arborescences
	  //des noeuds i et j à l'échelle courante+1.
	  VId root1 = _node1->getVertex();
	  TreeGraph* Tree1 = new TreeGraph(*(T1->getMTG()),root1,COMPO);
	  
	  VId root2 = _node2->getVertex();
	  TreeGraph* Tree2 = new TreeGraph(*(T2->getMTG()),root2,COMPO);
	  
	  MS_O_Matching* M=new MS_O_Matching(*Tree1,*Tree2,*CostMatrix,_init);
	  M->match();
	  alignement->add(M->getSequence());
	  matchCost += M->getMatchingCost();
	  insertCost += M->getInsertionCost();
	  deleteCost += M->getDeletionCost();
	  alignement->putNbIns(alignement->getNbIns()+M->getSequence()->getNbIns());
	  alignement->putNbDel(alignement->getNbDel()+M->getSequence()->getNbDel()); 
	  delete (MS_O_Matching*) M;
	  delete (TreeGraph*) Tree1;
	  delete (TreeGraph*) Tree2;
	  
	}
      subAlign(li,i-1,lj,j-1);
      
      return;	
    }
  }
  else
    {
      //sinon foret[li...i] et foret[lj...j] ne sont pas des arbres      
      DistanceType newmatch = fDist[fNum][idx(li,L1[i]-1)][idx(lj,L2[j]-1)];     
      
      if (f - newmatch - getTreeDistance(i,j) == 0) {
	subAlign(L1[i],i,L2[j],j);
	subAlign(li,L1[i]-1,lj,L2[j]-1);
	return;	
      }
      
      if (f - fi_1 - delLocalNode(T1,i) == 0) {
	subAlign(li,i-1,lj,j);
	deleteCost+=delLocalNode(T1,i);
      alignement->putNbDel(alignement->getNbDel()+getNbVertexAtScaleInit(T1, i)); 
	return ;
      } 
      
      if (f - fj_1 - insLocalNode(T2,j) == 0) {
	subAlign(li,i,lj,j-1);
	insertCost +=insLocalNode(T2,j);    
	alignement->putNbIns(alignement->getNbIns()+getNbVertexAtScaleInit(T2, j)); 
	return;
      } 
    }
}

int MS_O_Matching::getKey1(int li){
  int i;
  for (i = sizeKeyroots1 - 1; L1[keyroots1[i]] != L1[li]; i--);
  return keyroots1[i];
  
}

int MS_O_Matching::getKey2(int lj){
  int j;
  for (j = sizeKeyroots2 - 1; L2[keyroots2[j]] != L2[lj]; j--);
  return keyroots2[j];
}

int MS_O_Matching::getNbVertexAtScaleInit(TreeGraph* T, int i) {
  int nb = 0;
  int iReal = T->getRealNumber(i);
  TreeNode* _node =  T->getNode(iReal);
  if(T->getMTG()->vscale(_node->getVertex()) == _init){
    nb++;
  }
  else
    {
      VId root = _node->getVertex();
      TreeGraph* Tree = new TreeGraph(*(T->getMTG()),root,COMPO);
      int nbVertex = Tree->getNbVertex();
      for(int k = 0; k < nbVertex;k++)
	nb += getNbVertexAtScaleInit(Tree, k);
    }
  return nb;
}

