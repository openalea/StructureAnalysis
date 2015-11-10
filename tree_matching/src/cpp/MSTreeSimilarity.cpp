#include "MSSimilarity.h"
#include <iostream>
#include <cstdio>
#include <cmath>




// -------------
// Constructeur
// -------------
MSSimilarity::MSSimilarity(TreeGraph& input,TreeGraph& reference, NodeCost& nodeDistance, int init)
{
  T1=&input;
  T2=&reference;
  
  int sizeT1 = T1->getNbVertex();
  int sizeT2 = T2->getNbVertex();
  
  //Fonction de coût  
  CostMatrix = &nodeDistance;

  //flag indiquant le premier objet MSSimilarity
  _init = init;

  //tableau de distances entre arborescences
  _dst1.resize(sizeT1);
  int i= 0;
  for( i=0;i<=sizeT1-1;i++)
    _dst1[i].resize(sizeT2);
  
  //tableau de distance entre forêts
   int sizeFD = sizeT1*sizeT2;   
   fDist.resize(sizeFD);
   for( i = 0; i < sizeFD; i++){
     fDist[i].resize(1+sizeT1);
     for( int k = 0; k < 1+sizeT1; k++){
       fDist[i][k].resize(sizeT2+1);
     }
   }
   
   //sequence correspondant à l'alignement
   alignement = new Sequence(); 

   //Coûts associés à l'alignement trouvé
   insertCost = 0;
   deleteCost = 0;
   matchCost = 0; 
   
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
MSSimilarity::~MSSimilarity()
{
   
}

const float MSSimilarity::epsilon = 0.01;

float MSSimilarity::CostMatrix1[NB_SYMBOLS][NB_SYMBOLS]
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





// ----------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
void MSSimilarity::computeTreeDist(int input_vertex,int reference_vertex)
{

  //Calcul récursive des distances entre les sous-forêts
  computeForestsDistances(input_vertex, reference_vertex);
  float d = fDist[fix(input_vertex, reference_vertex)][idx(L1[input_vertex],input_vertex)][idx(L2[reference_vertex],reference_vertex)]; 
 _dst1[input_vertex][reference_vertex]=d;
}

//------------------------------------------------------------------------------------
//Distance entre deux arbres indexés par leurs racines
//-------------------------------------------------------------------------------------
float MSSimilarity::getTreeDistance(int input_vertex,int reference_vertex)
{
  
    return _dst1[input_vertex][reference_vertex];
}



//---------------------------------------------------------------------------------------
// Calcul de la distance entre les deux arbres T1 et T2
//---------------------------------------------------------------------------------------
DistanceType  MSSimilarity::match()
{
  int i,j,v,w;
  float dist;
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
      //  if(T1->getMTG()->vscale(T1->getNode(v)->getVertex()) < _init){
      //cout<<v<<" "<<w<<" "<<getTreeDistance(v,w)<<endl;}
    }
    
    
  }
  //dist = _dst1[keyroots1[sizeKeyroots1-1]][keyroots2[sizeKeyroots2-1]];
  
  //Calcul d'un alignement optimal
  makeAlignement(keyroots1[size1-1],keyroots2[size2-1]);
  dist = insertCost + deleteCost + matchCost;
  
  return dist;
}

void MSSimilarity::printForest(int v, int w){
  int fNum = fix(v,w);
  int i;
  for (i = L1[v]; i <= v; i++){
    cout<<"\t"<<i;
  }
  cout<<endl;
  int j;
  for (j = L2[w]; j <= w; j++){
    cout<<j;
    for (i = L1[v]; i <= v; i++){
      cout<<"\t"<<fDist[fNum][idx(L1[v],i)][idx(L2[w],j)];
    }
    cout<<endl;
  }
}

//-----------------------------------------------------------
// Calcul d'un alignement optimal entre les arborescences de 
// racines étiquetés i et j suivant la numérotation postfixe
//------------------------------------------------------------ 

void MSSimilarity::makeAlignement(int i, int j){
  
  int li = L1[i];
  int lj = L2[j];
  subAlign(li, i, lj, j);
}

//-----------------------------------------------------------
// Calcul d'un alignement optimal entre les forëts de 
// T1[li..i] et T2[lj..j] suivant la numérotation postfixe
//------------------------------------------------------------ 

void MSSimilarity::subAlign(int li, int i, int lj, int j){
  int num10;
  int num11;
  int  num12;
  int  num20;
  int num21;
  int num22;
  
 int num1,num2,num3,num4;
  
  if (L1[i] < li || L2[j] < lj) {
    return;
  }
  
  if(i==-1){
    for(int k=lj; k<=j; k++){
      insertCost +=insLocalNode(T2,k);
    }
    return;
  }
  
  if(j==-1){
   for(int k=li; k<=i; k++){
     deleteCost +=delLocalNode(T1,k);
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
      
      alignement->append( _node1->getVertex(),_node2->getVertex(), matchLocalNode(i,j));   
      /*const Feature* f1 = T1->getMTG()->si_feature(_node1->getVertex(),"base1id");
      const Feature* f2 = T2->getMTG()->si_feature(_node2->getVertex(),"base1id");
      alignement->append( f1->i,f2->i, matchLocalNode(i,j)); 
      const Feature* f3 = T1->getMTG()->si_feature(_node1->getVertex(),"base2id");
      const Feature* f4 = T2->getMTG()->si_feature(_node2->getVertex(),"base2id");
      alignement->append( f3->i,f4->i, matchLocalNode(i,j));*/ 
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
	
	MSSimilarity* M=new MSSimilarity(*Tree1,*Tree2,*CostMatrix,_init);
	M->match();
	alignement->add(M->getSequence());
	matchCost += M->getMatchingCost();
	insertCost += M->getInsertionCost();
	deleteCost += M->getDeletionCost();
	
	delete (MSSimilarity*) M;
	delete (TreeGraph*) Tree1;
	delete (TreeGraph*) Tree2;
      }
    return;
  } 
  
  
  int fNum = fix(getKey1(li),getKey2(lj));
  float f =  fDist[fNum][idx(li,i)][idx(lj,j)];
  float fi_1j_1 = fDist[fNum][idx(li,i-1)][idx(lj,j-1)];
  float fi_1 = fDist[fNum][idx(li,i-1)][idx(lj,j)];
  float fj_1 = fDist[fNum][idx(li,i)][idx(lj,j-1)];
  
  
  if(L1[i]==li && L2[j]==lj) {
    // si foret[li...i] et foret[lj...j] sont des arbres 
    
    if (f - fi_1 - delLocalNode(T1,i) == 0) {
      subAlign(li,i-1,lj,j);
      deleteCost+=delLocalNode(T1,i);      
      return ;
    } 
    
    if (f - fj_1 - insLocalNode(T2,j) == 0) {
      subAlign(li,i,lj,j-1);
      insertCost +=insLocalNode(T2,j);
      return;
    } 
    
    if (f - fi_1j_1 - matchLocalNode(i,j) == 0) {
      
      if(T1->getMTG()->vscale(_node1->getVertex()) == _init){
	//Si il s'agit de la dernière échelle de comparaison, on ajoute le couple de noeuds
	
	alignement->append( _node1->getVertex(),_node2->getVertex(), matchLocalNode(i,j));
	/*const Feature* f1 = T1->getMTG()->si_feature(_node1->getVertex(),"base1id");
      const Feature* f2 = T2->getMTG()->si_feature(_node2->getVertex(),"base1id");
      alignement->append( f1->i,f2->i, matchLocalNode(i,j)); 
      const Feature* f3 = T1->getMTG()->si_feature(_node1->getVertex(),"base2id");
      const Feature* f4 = T2->getMTG()->si_feature(_node2->getVertex(),"base2id");
      alignement->append( f3->i,f4->i, matchLocalNode(i,j));*/ 
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
	  
	  MSSimilarity* M=new MSSimilarity(*Tree1,*Tree2,*CostMatrix,_init);
	  M->match();
	  alignement->add(M->getSequence());
	  matchCost += M->getMatchingCost();
	  insertCost += M->getInsertionCost();
	  deleteCost += M->getDeletionCost();
	  
	  delete (MSSimilarity*) M;
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
      float newmatch = fDist[fNum][idx(li,L1[i]-1)][idx(lj,L2[j]-1)];     
      
      if (f - newmatch - _dst1[i][j] == 0) {
	subAlign(L1[i],i,L2[j],j);
	subAlign(li,L1[i]-1,lj,L2[j]-1);
	return;	
      }
      
      if (f - fi_1 - delLocalNode(T1,i) == 0) {
	subAlign(li,i-1,lj,j);
	deleteCost+=delLocalNode(T1,i);
	return ;
      } 
      
      if (f - fj_1 - insLocalNode(T2,j) == 0) {
	subAlign(li,i,lj,j-1);
	insertCost +=insLocalNode(T2,j);    
	return;
      } 
    }
}

int MSSimilarity::getKey1(int li){
  int i;
  for (i = sizeKeyroots1 - 1; L1[keyroots1[i]] != li; i--);
  return keyroots1[i];
  
}

int MSSimilarity::getKey2(int lj){
  int j;
  for (j = sizeKeyroots2 - 1; L2[keyroots2[j]] != lj; j--);
  return keyroots2[j];
}
