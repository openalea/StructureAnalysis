/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/matching_GPOSDE.cpp,v $
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


#include "matching_GPOSDE.h"

const int EMPTY_CLASS=-1;


// -------------
  // Constructeur
  // -------------
MatchingGposde::MatchingGposde(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  int i;
  T1=&input;
  T2=&reference;
  CostMatrix=&nodeDistance;
  int Size1 = T1->getNbVertex();
  int Size2 = T2->getNbVertex();
  _delete.resize(Size1);
  _insert.resize(Size2);
  
  _distances = new MatchingDistanceTable(*T1,*T2,nodeDistance);
  
  Size1 = T1->getNbClass();
  Size2 = T2->getNbClass();
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
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  _choicesFClass.resize(T1->getNbClass(),T2->getNbClass());

}
// -------------

// Destructeur
  // -------------
MatchingGposde::~MatchingGposde()
{
  delete (Sequence*) _alignement;
  delete (MatchingDistanceTable*) _distances;
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
void MatchingGposde::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  int i;
  DistanceType min;
 
  if(input_vertex!=EMPTY_NODE && reference_vertex==EMPTY_NODE){
    min = getDeleteCost(input_vertex) + getDBFC(T1->minimalClass(input_vertex),EMPTY_CLASS);
    _delete[input_vertex] = min;
    _distances->putInputTreeToEmpty(input_vertex,min);
  }
  
  if(input_vertex==EMPTY_NODE && reference_vertex!=EMPTY_NODE){
    min = getInsertCost(reference_vertex) + getDBFC(EMPTY_CLASS,T2->minimalClass(reference_vertex));
    _insert[reference_vertex]=min;
    _distances->putReferenceTreeFromEmpty(reference_vertex,min);
  }
  
  if(input_vertex!=EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    {
      int im=0,jm=0,MTC=0;
      
      DistanceType dist;
      
      int classmini = T1->minimalClass(input_vertex);
      int classminj = T2->minimalClass(reference_vertex);
      
      
      
      DistanceType nodeinodej = getDBFC(classmini,classminj) + getSubstitutionCost( input_vertex,reference_vertex);
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
      
      _choices.createList(input_vertex,reference_vertex);

      switch (MTC)
	{
	case 3 :{
	  _choices.putFirst(input_vertex,reference_vertex,im);
	  insert=getInBT(im,reference_vertex);
	  del=getDBT(input_vertex,EMPTY_TREE)-getDBT(im,EMPTY_TREE)+getDeBT(im,reference_vertex);
	  sub=getSuBT(im,reference_vertex);
	}break;
	case 2 :{
	  _choices.putFirst(input_vertex,reference_vertex,jm);
	  insert=getDBT(EMPTY_TREE,reference_vertex)-getDBT(EMPTY_TREE,jm)+getInBT(input_vertex,jm);
	  del=getDeBT(input_vertex,jm);
	  sub=getSuBT(input_vertex,jm);
	}break;
	case 1 :{
	  _choices.putFirst(input_vertex,reference_vertex,-1);
	  insert=getInBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
	del=getDeBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
	sub=getSuBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex))+  getSubstitutionCost(input_vertex,reference_vertex);
      }break;
	default :   assert(0);break;
	}
      _choices.putFirst(input_vertex,reference_vertex,MTC);
      // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
      // reference_vertex.
      _insertCost->putDBT(input_vertex,reference_vertex,insert);
      _deleteCost->putDBT(input_vertex,reference_vertex,del);
      _substitutionCost->putDBT(input_vertex,reference_vertex,sub);
      
    }
}


void MatchingGposde::distanceBetweenForestClass(int input_class,int reference_class)
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
      int im=0,jm=0,MFC=0;
      int rightclassi = T1->rightClasse(input_class);
      int rightclassj = T2->rightClasse(reference_class);
      
      DistanceType classi_classj = getDBC(input_class ,reference_class) + getDBFC(rightclassi,rightclassj);
      min = classi_classj; MFC = 1;
      
      DistanceType  classj_emptyclass =  getDBFC(input_class,rightclassj) + getDBC(EMPTY_CLASS,reference_class);
      if(min > classj_emptyclass) {min = classj_emptyclass;MFC=2;}
      
      DistanceType  classi_emptyclass =  getDBFC(rightclassi,reference_class) + getDBC(input_class,EMPTY_CLASS);
      if(min > classi_emptyclass) {min = classi_emptyclass;MFC=3;}
      
      DistanceType  emptytree_forestj =getDBFC(EMPTY_CLASS,reference_class);
      NodeList* nj  = T2->getRootsInForestClass(reference_class);
      vector<int>::const_iterator begin_node;
      begin_node = nj->begin();
      while (begin_node !=  nj->end())
	{
	  int reference_root=*begin_node;
	  DistanceType dist=getDBFC(input_class,T2->minimalClass(reference_root))-getDBT(EMPTY_TREE,reference_root) + emptytree_forestj + getInsertCost(reference_root);
	  if (min > dist) {min=dist;MFC=4;jm=reference_root;}
	  begin_node++;
	}
      
      DistanceType  emptytree_foresti =getDBFC(input_class,EMPTY_CLASS);
      NodeList* ni  = T1->getRootsInForestClass(input_class);
      
      begin_node = ni->begin();
      while (begin_node !=  ni->end())
	{
	  int input_root=*begin_node;
	  DistanceType dist=getDBFC(T1->minimalClass(input_root),reference_class)-getDBT(input_root,EMPTY_TREE) + emptytree_foresti + getDeleteCost(input_root);
	  if (min > dist) {min=dist;MFC=5;im=input_root;}
	  begin_node++;
	}

      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      switch(MFC)
	{
	case 5 :
	  {
	    _choicesFClass.putFirst(input_class,reference_class,T1->minimalClass(im));
	    del=getDBFC(input_class,EMPTY_CLASS)-getDBFC(T1->minimalClass(im),EMPTY_CLASS)+getDeBFC(T1->minimalClass(im),reference_class);
	    insert=getInBFC(T1->minimalClass(im),reference_class);
	    sub=getSuBFC(T1->minimalClass(im),reference_class);
	  }
	  break;
	case 4 :
	  {
	    _choicesFClass.putFirst(input_class,reference_class,T2->minimalClass(jm));
	    del=getDeBFC(input_class,T2->minimalClass(jm));
	    insert=getDBFC(EMPTY_TREE,reference_class)-getDBFC(EMPTY_TREE,T2->minimalClass(jm))+getInBFC(input_class,T2->minimalClass(jm));
	    sub=getSuBFC(input_class,T2->minimalClass(jm));
	  }
	  break;
	case 3 :
	  {
	    _choicesFClass.putFirst(input_class,reference_class,T1->rightClasse(input_class));
	    del=getDBC(input_class,EMPTY_CLASS) + getDeBFC(T1->rightClasse(input_class),reference_class);
	    insert= getInBFC(T1->rightClasse(input_class),reference_class);
	    sub = getSuBFC(T1->rightClasse(input_class),reference_class);
	  }
	  break;
	case 2 :
	  {
	    _choicesFClass.putFirst(input_class,reference_class,T2->rightClasse(reference_class));
	    del =  getDeBFC(input_class,T2->rightClasse(reference_class));
	    
	    insert= getDBC(EMPTY_CLASS,reference_class) + getDeBFC(input_class,T2->rightClasse(reference_class));
	    sub = getSuBFC(input_class,T2->rightClasse(reference_class));
	  }
	  break;
	case 1 :
	  {
	    _choicesFClass.putFirst(input_class,reference_class,-1);
	    del =  getDeBC(T1->rightClasse(input_class),T2->rightClasse(reference_class)) + getDeBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	    
	    insert=  getInBC(T1->rightClasse(input_class),T2->rightClasse(reference_class)) + getInBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	    sub =  getSuBC(T1->rightClasse(input_class),T2->rightClasse(reference_class)) + getSuBFC(T1->rightClasse(input_class),T2->rightClasse(reference_class));
	  }
	  break;
	default :   break;
	}
       _choicesFClass.putFirst(input_class,reference_class,MFC);
      _insertClassCost->putDBF(input_class,reference_class,insert);
      _deleteClassCost->putDBF(input_class,reference_class,del);
      _substitutionClassCost->putDBF(input_class,reference_class,sub);
    }
  _forestclasses[input_class+1][reference_class+1] = min;
  
}



void MatchingGposde::distanceBetweenClass(int input_class,int reference_class){

  int i;
  DistanceType dist =0;
  vector<int>::const_iterator begin_node;
  int ni = T1->getNbRootsInClass(input_class);
  int nj = T2->getNbRootsInClass(reference_class);
  NodeList* input_list= T1->getRootsInClass(input_class);
  NodeList* reference_list=T2->getRootsInClass(reference_class);
  MatchPath _restrMapp;
  
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
	    _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*(_distances->getDistanceTable())); 
	    _restrMapp.make(*input_list,*reference_list);

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
	    
	    
	    dist=_restrMapp.minCostFlow(_restrMappList);  
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
	_choicesFClass.createList(input_class,reference_class);
	sub=0;
	del=0;
	int i=1;
	NodeList* n = T1->getRootsInClass(input_class);
	begin_node = n->begin();
	while (begin_node !=  n->end())
	  {
	    int input_root=*begin_node;
	    _choicesFClass.putLast(input_class,reference_class,_restrMapp.who(_restrMappList[i]));
	    
	    if (_restrMapp.who(_restrMappList[i])!=EMPTY_NODE)
	      {
		sub=sub+getSuBT(input_root,_restrMapp.who(_restrMappList[i]));
		del=del+getDeBT(input_root,_restrMapp.who(_restrMappList[i]));
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
}

//Operator
DistanceType MatchingGposde::getDBT(int input_vertex,int reference_vertex) {
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

DistanceType MatchingGposde::getDBFC(int input_class,int reference_class) {
  return _forestclasses[input_class+1][reference_class+1];
} 
DistanceType MatchingGposde::getDBC(int input_class,int reference_class) {
  return _classes[input_class+1][reference_class+1];
} 

DistanceType MatchingGposde::getInsertCost(int input_vertex){
  return 1;
} 
DistanceType MatchingGposde::getDeleteCost(int reference_vertex){
  return 1;
} 

DistanceType MatchingGposde::getSubstitutionCost(int input_vertex,int reference_vertex){
  return 0;
} 

DistanceType  MatchingGposde::match(){
  int i,j;
  int size1 = T1->getNbClass();
  int size2 = T2->getNbClass();
  vector<int>::const_iterator begin_node1;
  vector<int>::const_iterator begin_node2;
  
  
  
  distanceBetweenClass(EMPTY_CLASS,EMPTY_CLASS);
  distanceBetweenForestClass(EMPTY_CLASS,EMPTY_CLASS);
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

DistanceType MatchingGposde::getInBT(int input_vertex,int reference_vertex) const
{
  return(_insertCost->getDBT(input_vertex,reference_vertex));
}

DistanceType MatchingGposde::getInBC(int input_class,int reference_class) const
{
  return(_insertClassCost->getDBT(input_class,reference_class));
}

DistanceType MatchingGposde::getInBFC(int input_class,int reference_class) const
{
  return(_insertClassCost->getDBF(input_class,reference_class));
}

DistanceType MatchingGposde::getDeBT(int input_vertex,int reference_vertex) const
{
  return(_deleteCost->getDBT(input_vertex,reference_vertex));
}

DistanceType MatchingGposde::getDeBC(int input_class,int reference_class) const
{
  return(_deleteClassCost->getDBT(input_class,reference_class));
}


DistanceType MatchingGposde::getDeBFC(int input_class,int reference_class) const
{
  return(_deleteClassCost->getDBF(input_class,reference_class));
}

DistanceType MatchingGposde::getSuBT(int input_vertex,int reference_vertex) const
{
  return(_substitutionCost->getDBT(input_vertex,reference_vertex));
}

DistanceType MatchingGposde::getSuBC(int input_class,int reference_class) const
{
  return(_substitutionClassCost->getDBT(input_class,reference_class));
}


DistanceType MatchingGposde::getSuBFC(int input_class,int reference_class) const
{
  return(_substitutionClassCost->getDBF(input_class,reference_class));
}




// renvoie le dernier element de la liste de la case node du tableau maintenant les listes d'alignement
int MatchingGposde::M(int i_node,int r_node)
{
  return(_choices.getList(i_node,r_node)->back());
}

int MatchingGposde::Lat(ChoiceList* L, int vertex)
{
  ChoiceList::iterator begin;
  begin = L->begin();
  for (int i=0;i<vertex;i++)
    begin++;
  return(*begin);
}



void MatchingGposde::TreeList(int input_vertex,int reference_vertex,Sequence& sequence)
{
  if ((!T1->isNull())&&(!T2->isNull()))
    {
      if ((input_vertex!=-1)&&(reference_vertex!=-1))
	{
	  ChoiceList* L=_choices.getList(input_vertex,reference_vertex);
	  int tree_choice=L->front();
	  switch(tree_choice)
	    {
	    case 1:
	      {
		sequence.append(input_vertex,reference_vertex,getSubstitutionCost(input_vertex,reference_vertex));
		ForestClassList(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex),sequence);
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

void MatchingGposde::ForestClassList(int input_class,int reference_class,Sequence& sequence)
{
  if ((input_class!=-1)&&(reference_class!=-1))
    {
      ChoiceList* L=_choicesFClass.getList(input_class,reference_class);
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
	case 4: ForestClassList(input_class,Lat(L,1),sequence);break;
	case 5: ForestClassList(Lat(L,1),reference_class,sequence);break;
	default : break;
	}
    }
}
