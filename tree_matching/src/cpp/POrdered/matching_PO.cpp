/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/POrdered/matching_PO.cpp,v $
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


#include "matching_PO.h"

const int EMPTY_CLASS=-1;


// -------------
  // Constructeur
  // -------------
Matching_PO::Matching_PO(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
	int i;
  _i_v=0;
  _r_v=0;
  T1=&input;
  T2=&reference;
  CostMatrix=&nodeDistance;
  int Size1 = T1->getNbVertex();
  int Size2 = T2->getNbVertex();
  _delete.resize(Size1);
  _insert.resize(Size2);
  
  _distances = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  
  Size1 = T1->getNbClass();
  sizeC1=Size1;
  Size2 = T2->getNbClass();
  sizeC2=Size2;
  
  _forestclasses.resize(Size1+1);
  for(i =0; i <= Size1;i++)
    _forestclasses[i].resize(Size2+1);
  
  _classes.resize(Size1+1);
  for(i =0; i <= Size1;i++)
    _classes[i].resize(Size2+1);
  
  
  _alignement = new Sequence(); 
  
  _insertCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  _deleteCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  _substitutionCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  _insertClassCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  _deleteClassCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  _substitutionClassCost = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  
  // _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
  // c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices = new ChoiceTable(T1->getNbVertex(),T2->getNbVertex());
  _choicesFClass =  new ChoiceTable(T1->getNbClass(),T2->getNbClass());

}
// -------------

// Destructeur
  // -------------
Matching_PO::~Matching_PO()
{
	int i;
  delete (Sequence*) _alignement;
  delete (MatchingDistanceTable*) _distances;
   delete (MatchingDistanceTable*) _insertCost;
   delete (MatchingDistanceTable*) _deleteCost;
   delete (MatchingDistanceTable*) _substitutionCost;
   delete (MatchingDistanceTable*) _insertClassCost;
   delete (MatchingDistanceTable*) _deleteClassCost;
   delete (MatchingDistanceTable*) _substitutionClassCost;
delete (ChoiceTable*) _choices;
delete (ChoiceTable*) _choicesFClass;

   _delete.clear();
  _insert.clear();
  for(i =0; i <= sizeC1;i++)
    _forestclasses[i].clear();
_forestclasses.clear();
  
  for(i =0; i <= sizeC1;i++)
    _classes[i].clear();
  _classes.clear();
  
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
void Matching_PO::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  int i;
  DistanceType min;
 
  if(input_vertex!=EMPTY_NODE && reference_vertex==EMPTY_NODE){
    min = getDeleteCost(input_vertex) + getDBF(input_vertex,EMPTY_NODE);
    _delete[input_vertex] = min;
    _distances->putInputTreeToEmpty(input_vertex,min);
  }
  
  if(input_vertex==EMPTY_NODE && reference_vertex!=EMPTY_NODE){
    min = getInsertCost(reference_vertex) + getDBF(EMPTY_NODE,reference_vertex);
    _insert[reference_vertex]=min;
    _distances->putReferenceTreeFromEmpty(reference_vertex,min);
  }
  
  if(input_vertex!=EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    {
      int im=0,jm=0,MTC=0;
      
      DistanceType dist;
      
      
      DistanceType nodeinodej = getDBF( input_vertex,reference_vertex) + getSubstitutionCost( input_vertex,reference_vertex);
      min=nodeinodej;
      //cout<<"nodeinodej "<<nodeinodej<<endl;
      MTC = 1;
      DistanceType  emptytree_treej =getDBT(EMPTY_TREE,reference_vertex);
      int nj = T2->getNbChild(reference_vertex);
      for (i=1;i<=nj;i++)
	{
	  int reference_child=T2->child(reference_vertex,i);
	  dist=getDBT(input_vertex,reference_child)-getDBT(EMPTY_TREE,reference_child) + emptytree_treej;
	  if (min > dist) {min=dist; MTC=2; jm =reference_child; }
	}
      //cout<<"emptytree_treej "<<dist<<endl;

      DistanceType  emptytree_treei =getDBT(input_vertex,EMPTY_TREE);
      int ni = T1->getNbChild(input_vertex);
      for (i=1;i<=ni;i++)
	{
	  int input_child=T1->child(input_vertex,i);
	  dist=getDBT(input_child,reference_vertex)-getDBT(input_child, EMPTY_TREE) + emptytree_treei;
	  if (min > dist) {min=dist; MTC=3; im = input_child;}
	}
      //cout<<"emptytree_treei "<<dist<<endl;
      _distances->putDBT(input_vertex,reference_vertex,min);
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      
     

      switch (MTC)
	{
	case 3 :{
	  _choices->putFirst(input_vertex,reference_vertex,im);
	  insert=getInBT(im,reference_vertex);
	  del=getDBT(input_vertex,EMPTY_TREE)-getDBT(im,EMPTY_TREE)+getDeBT(im,reference_vertex);
	  sub=getSuBT(im,reference_vertex);
	}break;
	case 2 :{
	  _choices->putFirst(input_vertex,reference_vertex,jm);
	  insert=getDBT(EMPTY_TREE,reference_vertex)-getDBT(EMPTY_TREE,jm)+getInBT(input_vertex,jm);
	  del=getDeBT(input_vertex,jm);
	  sub=getSuBT(input_vertex,jm);
	}break;
	case 1 :{
	  _choices->putFirst(input_vertex,reference_vertex,-1);
	  insert=getInBF(input_vertex,reference_vertex);
	del=getDeBF(input_vertex,reference_vertex);
	sub=getSuBF(input_vertex,reference_vertex)+  getSubstitutionCost(input_vertex,reference_vertex);
      }break;
	default :   assert(0);break;
	}
      _choices->putFirst(input_vertex,reference_vertex,MTC);
      // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
      // reference_vertex.
      _insertCost->putDBT(input_vertex,reference_vertex,insert);
      _deleteCost->putDBT(input_vertex,reference_vertex,del);
      _substitutionCost->putDBT(input_vertex,reference_vertex,sub);
      
    }
}


void Matching_PO::distanceBetweenForest(int input_vertex,int reference_vertex)
{  
  int i;
  
  DistanceType min;
  if(input_vertex==EMPTY_NODE && reference_vertex==EMPTY_NODE)
    min = 0;
  if(input_vertex!=EMPTY_NODE && reference_vertex==EMPTY_NODE)
    min = getDBFC(T1->minimalClass(input_vertex),EMPTY_CLASS);
  
  if(input_vertex==EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    min = getDBFC(EMPTY_CLASS,T2->minimalClass(reference_vertex));
  
  if(input_vertex!=EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    {
      int im=0,jm=0,MFC=0;
      
      DistanceType classi_classj = getDBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
      min = classi_classj; MFC = 1;
      
      DistanceType  emptytree_forestj =getDBF(EMPTY_NODE,reference_vertex);
      int nj  = T2->getNbChild(reference_vertex);
      for (i=1;i<=nj;i++)
	{
	  int reference_root=T2->child(reference_vertex,i);
	  DistanceType dist=getDBF(input_vertex,reference_root)-getDBT(EMPTY_TREE,reference_root) + emptytree_forestj + getInsertCost(reference_root);
	  if (min > dist) {min=dist;MFC=4;jm=reference_root;}
	}
      
      DistanceType  emptytree_foresti =getDBF(input_vertex,EMPTY_NODE);
      int ni  = T1->getNbChild(input_vertex);
      
      for (i=1;i<=ni;i++)
	{
	  int input_root=T1->child(input_vertex,i);
	  DistanceType dist=getDBF(input_root,reference_vertex)-getDBT(input_root,EMPTY_TREE) + emptytree_foresti + getDeleteCost(input_root);
	  if (min > dist) {min=dist;MFC=5;im=input_root;}
	}

      _choices->createList(input_vertex,reference_vertex);
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      switch(MFC)
	{
	case 5 :
	  {
	    _choices->putFirst(input_vertex,reference_vertex,im);
	    del=getDBF(input_vertex,EMPTY_NODE)-getDBF(im,EMPTY_NODE)+getDeBF(im,reference_vertex);
	    insert=getInBF(im,reference_vertex);
	    sub=getSuBF(im,reference_vertex);
	  }
	  break;
	case 4 :
	  {
	    _choices->putFirst(input_vertex,reference_vertex,jm);
	    del=getDeBF(input_vertex,jm);
	    insert=getDBF(EMPTY_NODE,reference_vertex)-getDBF(EMPTY_TREE,jm)+getInBF(input_vertex,jm);
	    sub=getSuBF(input_vertex,jm);
	  }
	  break;
	case 1 :
	  {
	    _choices->putFirst(input_vertex,reference_vertex,-1);
	    del =  getDeBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
	    
	    insert=  getInBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
	    sub =  getSuBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
	  }
	  break;
	default :   break;
	}
       _choices->putFirst(input_vertex,reference_vertex,MFC);
      _insertCost->putDBF(input_vertex,reference_vertex,insert);
      _deleteCost->putDBF(input_vertex,reference_vertex,del);
      _substitutionCost->putDBF(input_vertex,reference_vertex,sub);


      _distances->putDBF(input_vertex,reference_vertex,min);
    }
  
 

  
}

void Matching_PO::distanceBetweenForestClass(int input_class,int reference_class)
{  
  int i;
  
  DistanceType min;
  if(input_class==EMPTY_CLASS && reference_class==EMPTY_CLASS)
    min = 0;
  if(input_class!=EMPTY_CLASS && reference_class==EMPTY_CLASS)
    min = getDBC(input_class,EMPTY_CLASS) + getDBFC(T1->rightClasse(input_class),EMPTY_CLASS);
  
  if(input_class==EMPTY_CLASS && reference_class!=EMPTY_CLASS)
    min = getDBC(EMPTY_CLASS,reference_class) + getDBFC(EMPTY_CLASS,T2->rightClasse(reference_class));
  
  if(input_class!=EMPTY_CLASS && reference_class!=EMPTY_CLASS)
    {
      int MFC=0;
      int rightclassi = T1->rightClasse(input_class);
      int rightclassj = T2->rightClasse(reference_class);
      
      DistanceType classi_classj = getDBC(input_class ,reference_class) + getDBFC(rightclassi,rightclassj);
      min = classi_classj; MFC = 1;
      
      DistanceType  classj_emptyclass =  getDBFC(input_class,rightclassj) + getDBC(EMPTY_CLASS,reference_class);
      if(min > classj_emptyclass) {min = classj_emptyclass;MFC=2;}
      
      DistanceType  classi_emptyclass =  getDBFC(rightclassi,reference_class) + getDBC(input_class,EMPTY_CLASS);
      if(min > classi_emptyclass) {min = classi_emptyclass;MFC=3;}
      
   
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      switch(MFC)
	{
	case 3 :
	  {
	    _choicesFClass->putFirst(input_class,reference_class,T1->rightClasse(input_class));
	    del=getDBC(input_class,EMPTY_CLASS) + getDeBFC(T1->rightClasse(input_class),reference_class);
	    insert= getInBFC(T1->rightClasse(input_class),reference_class);
	    sub = getSuBFC(T1->rightClasse(input_class),reference_class);
	  }
	  break;
	case 2 :
	  {
	    _choicesFClass->putFirst(input_class,reference_class,T2->rightClasse(reference_class));
	    del =  getDeBFC(input_class,T2->rightClasse(reference_class));
	    
	    insert= getDBC(EMPTY_CLASS,reference_class) + getInBFC(input_class,T2->rightClasse(reference_class));
	    sub = getSuBFC(input_class,T2->rightClasse(reference_class));
	  }
	  break;
	case 1 :
	  {
	    _choicesFClass->putFirst(input_class,reference_class,-1);
	    del =  getDeBC(input_class,reference_class) + getDeBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	    
	    insert=  getInBC(input_class,reference_class) + getInBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	    sub =  getSuBC(input_class,reference_class) + getSuBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	  }
	  break;
	default :   break;
	}
      _choicesFClass->putFirst(input_class,reference_class,MFC);
      _insertClassCost->putDBF(input_class,reference_class,insert);
      _deleteClassCost->putDBF(input_class,reference_class,del);
      _substitutionClassCost->putDBF(input_class,reference_class,sub);
    }
  _forestclasses[input_class+1][reference_class+1] = min;
  
}


void Matching_PO::distanceBetweenClass(int input_class,int reference_class){

  int i;
  DistanceType dist =0;
  vector<int>::const_iterator begin_node;
  int ni = T1->getNbRootsInClass(input_class);
  int nj = T2->getNbRootsInClass(reference_class);
  NodeList* input_list= T1->getRootsInClass(input_class);
  NodeList* reference_list=T2->getRootsInClass(reference_class);
  MatchPath* _restrMapp = new MatchPath();
  
  VertexVector _restrMappList; 
  _restrMappList.resize(ni+nj+3,EMPTY_NODE);  
  
  if(ni ==0 && nj ==0)
    dist=0;
  else{
    if (ni==0)
      {
	_restrMappList[1]=2;
	for (i=1;i<=nj;i++) { _restrMappList[i+1]=1; }
	NodeList* n = T2->getRootsInClass(reference_class);
	begin_node = n->begin();
	while (begin_node !=  n->end())
	  {
	    int reference_root=*begin_node;
	    //dist= dist+ getDBFC(EMPTY_CLASS,T2->minimalClass(reference_root)) + getInsertCost(reference_root);
	    dist= dist+ getDBT(EMPTY_NODE,reference_root);
	    begin_node++;
	  }
      }
    else
      {
	if (nj==0){
	  _restrMappList[2]=1;
	  for (i=1;i<=ni;i++) { _restrMappList[i]=ni+1; }
	  NodeList* n = T1->getRootsInClass(input_class);
	  begin_node = n->begin();
	  while (begin_node !=  n->end())
	    {
	      int input_root=*begin_node;
	      //dist= dist+ getDBFC(T1->minimalClass(input_root),EMPTY_CLASS) + getDeleteCost(input_root);
	      dist= dist+ getDBT(input_root,EMPTY_NODE);
	      begin_node++;
	    }
	}
	
	else
	  {
	    _restrMapp->link(I_MAX(T1->getDegree(),T2->getDegree()),*(_distances->getDistanceTable())); 
	    _restrMapp->make(*input_list,*reference_list);
	    vector<int>::const_iterator begin_node;
	    //cout<<"input_list : ";
	    //begin_node = input_list->begin();
	    //while(begin_node != input_list->end()){
	    //cout<<*begin_node;
	    //begin_node++;
	    //}
	    //cout<<endl;
	    //cout<<"reference_list : ";
	    //begin_node = reference_list->begin();
	    //while(begin_node != reference_list->end()){
	    //cout<<*begin_node;
	    //begin_node++;
	    //}
	    //cout<<endl;
	    
	    
	    dist=_restrMapp->minCostFlow(_restrMappList);  
	    //dist = 0;
	  }
      }
  } 

  DistanceType del,insert,sub;
  del=0;insert=0;sub=0;
  
  if(ni!=0 || nj !=0){
    if (ni==0)
      {
	insert=getDBC(EMPTY_CLASS,reference_class);
	del=0;
	sub=0;
      }
    else if (nj==0)
      {
	insert=0;
	del=getDBC(input_class,EMPTY_CLASS);
	sub=0;
      }
    else
      {
	_choicesFClass->createList(input_class,reference_class);
	sub=0;
	del=0;
	int i=1;
	NodeList* n = T1->getRootsInClass(input_class);
	begin_node = n->begin();
	while (begin_node !=  n->end())
	  {
	    int input_root=*begin_node;
	    _choicesFClass->putLast(input_class,reference_class,_restrMapp->who(_restrMappList[i]));
	    
	    if (_restrMapp->who(_restrMappList[i])!=EMPTY_NODE)
	      {
		sub=sub+getSuBT(input_root,_restrMapp->who(_restrMappList[i]));
		del=del+getDeBT(input_root,_restrMapp->who(_restrMappList[i]));
	      }
	    else
	      {
		del=del+getDBT(input_root,EMPTY_TREE);
	      }
	    i++;
	    begin_node++;
	  }
	
	insert=dist-sub-del;
      } 
    if(input_class!=-1 && reference_class !=-1){
      _insertClassCost->putDBT(input_class,reference_class,insert);
      _deleteClassCost->putDBT(input_class,reference_class,del);
      _substitutionClassCost->putDBT(input_class,reference_class,sub);
    }
  }
  _classes[input_class+1][reference_class+1] = dist;
_restrMappList.clear();
delete (MatchPath*) _restrMapp;
}

//Operator
DistanceType Matching_PO::getDBT(int input_vertex,int reference_vertex) {
  DistanceType dist;
  if(input_vertex == EMPTY_NODE &&reference_vertex== EMPTY_NODE)
    dist = 0;
  if(input_vertex != EMPTY_NODE &&reference_vertex== EMPTY_NODE)
    dist = _delete[input_vertex];
  if(input_vertex == EMPTY_NODE &&reference_vertex!= EMPTY_NODE)
    dist = _insert[reference_vertex];
  if(input_vertex != EMPTY_NODE &&reference_vertex!= EMPTY_NODE)
    dist = _distances->getDBT(input_vertex,reference_vertex);  
  return dist;
} 

DistanceType Matching_PO::getDBF(int input_vertex,int reference_vertex) {
  DistanceType dist;
  if(input_vertex == EMPTY_NODE &&reference_vertex== EMPTY_NODE)
    dist = 0;
  if(input_vertex != EMPTY_NODE &&reference_vertex== EMPTY_NODE)
    dist = getDBFC(T1->minimalClass(input_vertex),EMPTY_CLASS);
  if(input_vertex == EMPTY_NODE && reference_vertex!= EMPTY_NODE)
    dist = getDBFC(EMPTY_CLASS,T2->minimalClass(reference_vertex));
  if(input_vertex != EMPTY_NODE &&reference_vertex!= EMPTY_NODE)
    dist = _distances->getDBF(input_vertex,reference_vertex);  
  return dist;
}

DistanceType Matching_PO::getDBFC(int input_class,int reference_class) {
  return _forestclasses[input_class+1][reference_class+1];
} 
DistanceType Matching_PO::getDBC(int input_class,int reference_class) {
  return _classes[input_class+1][reference_class+1];
} 

DistanceType Matching_PO::getInsertCost(int input_vertex){
  return _distances->getICost(input_vertex);
} 

DistanceType Matching_PO::getDeleteCost(int reference_vertex){
  return _distances->getDCost(reference_vertex);
} 

DistanceType Matching_PO::getSubstitutionCost(int input_vertex,int reference_vertex){
  return _distances->getCCost(input_vertex,reference_vertex);
} 

DistanceType  Matching_PO::match(){
 
  int i,j;
  int size1 = T1->getNbClass();
  int size2 = T2->getNbClass();
  vector<int>::const_iterator begin_node1;
  vector<int>::const_iterator begin_node2;
  
  
  
  distanceBetweenClass(EMPTY_CLASS,EMPTY_CLASS);
  distanceBetweenForestClass(EMPTY_CLASS,EMPTY_CLASS);
  distanceBetweenForest(EMPTY_NODE,EMPTY_NODE);
  distanceBetweenTree(EMPTY_NODE,EMPTY_NODE);
  
  //cout<<"delete"<<endl;
  for(i=size1-1;i>=0;i--){
    _insertClassCost->openDistancesVector(i);
    
    _deleteClassCost->openDistancesVector(i);
    
    _substitutionClassCost->openDistancesVector(i);
    
    NodeList* input_list= T1->getRootsInClass(i);
    begin_node1 = input_list->begin();
    while (begin_node1 !=  input_list->end())
      {
	int input_root=*begin_node1;
	//cout<<input_root<<endl;
	distanceBetweenForest(input_root,EMPTY_NODE);
	distanceBetweenTree(input_root,EMPTY_NODE);
	//cout<<"tree "<<input_root<<" "<<getDBT(input_root,EMPTY_NODE)<<endl;
	_distances->openDistancesVector(input_root);
	_insertCost->openDistancesVector(input_root);
	
	_deleteCost->openDistancesVector(input_root);
	
	_substitutionCost->openDistancesVector(input_root);
	begin_node1++;
      }
    distanceBetweenClass(i,EMPTY_CLASS);
    distanceBetweenForestClass(i,EMPTY_CLASS);
  }
  //cout<<"insert"<<endl;
  _distances->openDistancesVector(T1->getNbVertex());
  for(j=size2-1;j>=0;j--){
    NodeList* reference_list=T2->getRootsInClass(j);
    begin_node2 = reference_list->begin();
    while (begin_node2 !=  reference_list->end())
      {
	int reference_root=*begin_node2;
	//cout<<reference_root<<endl;
	distanceBetweenForest(EMPTY_NODE,reference_root);
	distanceBetweenTree(EMPTY_NODE,reference_root);
	//cout<<"tree "<<reference_root<<" "<<getDBT(EMPTY_NODE,reference_root)<<endl;
	begin_node2++;
      }
    distanceBetweenClass(EMPTY_CLASS,j);
    distanceBetweenForestClass(EMPTY_CLASS,j);
  } 

  //cout<<"match"<<endl;
  for(i=size1-1;i>=0;i--){
    for(j=size2-1;j>=0;j--){
      NodeList* input_list= T1->getRootsInClass(i);
      NodeList* reference_list=T2->getRootsInClass(j);
      begin_node1 = input_list->begin();
      while (begin_node1 !=  input_list->end())
	{
	  int input_root=*begin_node1;
	  begin_node2 = reference_list->begin();
	  while (begin_node2 !=  reference_list->end())
	    {
	      int reference_root=*begin_node2;
	      //cout<<input_root<<" ; "<<reference_root<<endl;
	      distanceBetweenForest(input_root,reference_root);
	      distanceBetweenTree(input_root,reference_root);
	      //cout<<"tree "<<input_root<<";"<<reference_root<<" "<<getDBT(input_root,reference_root)<<endl;
	      begin_node2++;
	    }
	  begin_node1++;
	}
      
      distanceBetweenClass(i,j);
      //cout<<"class "<<i<<";"<<j<<" "<<getDBC(i,j)<<endl;
      distanceBetweenForestClass(i,j);
      //cout<<"forestclass "<<i<<";"<<j<<" "<<getDBFC(i,j)<<endl;
    }
  }  
  _alignement->append(0,0,getSubstitutionCost(0,0));
  return getDBT(0,0);
} 

DistanceType Matching_PO::getInBT(int input_vertex,int reference_vertex) const
{
  return(_insertCost->getDBT(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getInBF(int input_vertex,int reference_vertex) const
{
  return(_insertCost->getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getInBC(int input_class,int reference_class) const
{
  return(_insertClassCost->getDBT(input_class,reference_class));
}

DistanceType Matching_PO::getInBFC(int input_class,int reference_class) const
{
  return(_insertClassCost->getDBF(input_class,reference_class));
}

DistanceType Matching_PO::getDeBT(int input_vertex,int reference_vertex) const
{
  return(_deleteCost->getDBT(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getDeBF(int input_vertex,int reference_vertex) const
{
  return(_deleteCost->getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getDeBC(int input_class,int reference_class) const
{
  return(_deleteClassCost->getDBT(input_class,reference_class));
}


DistanceType Matching_PO::getDeBFC(int input_class,int reference_class) const
{
  return(_deleteClassCost->getDBF(input_class,reference_class));
}

DistanceType Matching_PO::getSuBT(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost->getDBT(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getSuBF(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost->getDBF(input_vertex,reference_vertex));
}

DistanceType Matching_PO::getSuBC(int input_class,int reference_class) const
{
  return(_substitutionClassCost->getDBT(input_class,reference_class));
}


DistanceType Matching_PO::getSuBFC(int input_class,int reference_class) const
{
  return(_substitutionClassCost->getDBF(input_class,reference_class));
}




// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int Matching_PO::M(int i_node,int r_node)
{
  return(_choices->getList(i_node,r_node)->back());
}

int Matching_PO::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void Matching_PO::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      if ((input_vertex!=-1)&&(reference_vertex!=-1))
	{
	  ChoiceList* L=_choices->getList(input_vertex,reference_vertex);
	  int tree_choice=L->front();
	  switch(tree_choice)
	    {
	    case 1:
	      {
		sequence.append(input_vertex,reference_vertex,getSubstitutionCost(input_vertex,reference_vertex));
		ForestList(input_vertex,reference_vertex,sequence);
	      }
	      break;
	    case 2:
	      {
		TreeList(input_vertex,Lat(L,1),sequence);
	      }
	      break;
	    case 3: 
	      {	  
		TreeList(Lat(L,1),reference_vertex,sequence);
	      }
	      break;
	    default : break;
	    }
	}
    }
  
}

void Matching_PO::ForestList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      if ((input_vertex!=-1)&&(reference_vertex!=-1))
	{
	  ChoiceList* L=_choices->getList(input_vertex,reference_vertex);
	  int tree_choice=Lat(L,2);
	  switch(tree_choice)
	    {
	    case 1:
	      {
		ForestClassList(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex),sequence);
	      }
	      break;
	    case 4:
	      {
		ForestList(input_vertex,Lat(L,3),sequence);
	      }
	      break;
	    case 5: 
	      {	  
		ForestList(Lat(L,3),reference_vertex,sequence);
	      }
	      break;
	    default : break;
	    }
	}
    }
  
}


void Matching_PO::ForestClassList(int input_class,int reference_class,Sequence& sequence)
{
  if ((input_class!=-1)&&(reference_class!=-1))
    {
      ChoiceList* L=_choicesFClass->getList(input_class,reference_class);
      int forest_choice=L->front();
      switch(forest_choice)
	{
	case 1: 
	  {
	    int i =1;
	    NodeList* n = T1->getRootsInClass(input_class);
	    vector<int>::const_iterator begin_node;
	    begin_node = n->begin();
	    while (begin_node !=  n->end()){
	      int input_root=*begin_node;
	      int reference_root=Lat(L,1+i);
	      if (reference_root!=-1) 
		TreeList(input_root,reference_root,sequence);
	      i++;
	      begin_node++;
	    }
	    
	    ForestClassList(T1->rightClasse(input_class),T2->rightClasse(reference_class),sequence);
	  }
	  break;
	case 2: ForestClassList(input_class,Lat(L,1),sequence);break;
	case 3: ForestClassList(Lat(L,1),reference_class,sequence);break;
	default : break;
	}
    }
}
