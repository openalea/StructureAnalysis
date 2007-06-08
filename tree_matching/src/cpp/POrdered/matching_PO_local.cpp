/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (ouangrao@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/POrdered/matching_PO_local.cpp,v $
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


#include "matching_PO_local.h"

const int EMPTY_CLASS=-1;


// -------------
  // Constructeur
  // -------------
Matching_PO_Local::Matching_PO_Local(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  Matching_PO(input,reference,nodeDistance)
{
}

// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[input_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
void Matching_PO_Local::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  int i;
  DistanceType max;
 
  if(input_vertex!=EMPTY_NODE && reference_vertex==EMPTY_NODE){
    max = 0;
    _delete[input_vertex] = max;
    _distances->putInputTreeToEmpty(input_vertex,max);
  }
  
  if(input_vertex==EMPTY_NODE && reference_vertex!=EMPTY_NODE){
    max = 0;
    _insert[reference_vertex]=max;
    _distances->putReferenceTreeFromEmpty(reference_vertex,max);
  }
  
  if(input_vertex!=EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    {
      int im=0,jm=0,MTC=0;
      
      DistanceType dist;
      
      
      DistanceType nodeinodej = getDBF( input_vertex,reference_vertex) + getSubstitutionCost( input_vertex,reference_vertex);
      max=nodeinodej;
      //cout<<"nodeinodej "<<nodeinodej<<endl;
      MTC = 1;
      
      int nj = T2->getNbChild(reference_vertex);
      for (i=1;i<=nj;i++)
	{
	  int reference_child=T2->child(reference_vertex,i);
	  dist=getDBT(input_vertex,reference_child) + getInsertCost(reference_vertex);
	  if (max < dist) {max=dist; MTC=2; jm =reference_child; }
	}
      //cout<<"emptytree_treej "<<dist<<endl;

      int ni = T1->getNbChild(input_vertex);
      for (i=1;i<=ni;i++)
	{
	  int input_child=T1->child(input_vertex,i);
	  dist=getDBT(input_child,reference_vertex) + getDeleteCost(input_vertex);
	  if (max < dist) {max=dist; MTC=3; im = input_child;}
	}
      //cout<<"emptytree_treei "<<dist<<endl;
     
      _distances->putDBT(input_vertex,reference_vertex,max);
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      
     

      switch (MTC)
	{
	case 3 :{
	  _choices->putFirst(input_vertex,reference_vertex,im);
	  insert=getInBT(im,reference_vertex);
	  del=getDeleteCost(input_vertex)+getDeBT(im,reference_vertex);
	  sub=getSuBT(im,reference_vertex);
	}break;
	case 2 :{
	  _choices->putFirst(input_vertex,reference_vertex,jm);
	  insert=getInsertCost(reference_vertex)+getInBT(input_vertex,jm);
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


void Matching_PO_Local::distanceBetweenForest(int input_vertex,int reference_vertex)
{  
  int i;
  
  DistanceType   max;
  if(input_vertex==EMPTY_NODE && reference_vertex==EMPTY_NODE)
    max = 0;
  if(input_vertex!=EMPTY_NODE && reference_vertex==EMPTY_NODE)
    max = 0;
  
  if(input_vertex==EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    max = 0;
  
  if(input_vertex!=EMPTY_NODE && reference_vertex!=EMPTY_NODE)
    {
      int im=0,jm=0,MFC=0;
      
      DistanceType classi_classj = getDBFC(T1->minimalClass(input_vertex),T2->minimalClass(reference_vertex));
      max = classi_classj; MFC = 1;
      
      int nj  = T2->getNbChild(reference_vertex);
      for (i=1;i<=nj;i++)
	{
	  int reference_root=T2->child(reference_vertex,i);
	  DistanceType dist=getDBF(input_vertex,reference_root)+ getInsertCost(reference_root);
	  if (max < dist) {max=dist;MFC=4;jm=reference_root;}
	}
      

      int ni  = T1->getNbChild(input_vertex);
      
      for (i=1;i<=ni;i++)
	{
	  int input_root=T1->child(input_vertex,i);
	  DistanceType dist=getDBF(input_root,reference_vertex)+ getDeleteCost(input_root);
	  if (max < dist) {max=dist;MFC=5;im=input_root;}
	}

      _choices->createList(input_vertex,reference_vertex);
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      switch(MFC)
	{
	case 5 :
	  {
	    _choices->putFirst(input_vertex,reference_vertex,im);
	    del=getDeleteCost(input_vertex)+getDeBF(im,reference_vertex);
	    insert=getInBF(im,reference_vertex);
	    sub=getSuBF(im,reference_vertex);
	  }
	  break;
	case 4 :
	  {
	    _choices->putFirst(input_vertex,reference_vertex,jm);
	    del=getDeBF(input_vertex,jm);
	    insert=getInsertCost(reference_vertex)+getInBF(input_vertex,jm);
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


      _distances->putDBF(input_vertex,reference_vertex,max);
    }
  
 

  
}

void Matching_PO_Local::distanceBetweenForestClass(int input_class,int reference_class)
{  
  int i;
  
  DistanceType max;
  if(input_class==EMPTY_CLASS && reference_class==EMPTY_CLASS)
    max = 0;
  if(input_class!=EMPTY_CLASS && reference_class==EMPTY_CLASS)
    max = 0;
  
  if(input_class==EMPTY_CLASS && reference_class!=EMPTY_CLASS)
    max = 0;
  
  if(input_class!=EMPTY_CLASS && reference_class!=EMPTY_CLASS)
    {
      int MFC=0;
      int rightclassi = T1->rightClasse(input_class);
      int rightclassj = T2->rightClasse(reference_class);
      
      DistanceType classi_classj = getDBC(input_class ,reference_class) + getDBFC(rightclassi,rightclassj);
      max = classi_classj; MFC = 1;
      
      DistanceType  classj_emptyclass =  getDBFC(input_class,rightclassj);
      if(max < classj_emptyclass) {max = classj_emptyclass;MFC=2;}
      
      DistanceType  classi_emptyclass =  getDBFC(rightclassi,reference_class);
      if(max < classi_emptyclass) {max = classi_emptyclass;MFC=3;}
      
   
      DistanceType del,insert,sub;
      del=0;insert=0;sub=0;
      switch(MFC)
	{
	case 3 :
	  {
	    _choicesFClass->putFirst(input_class,reference_class,T1->rightClasse(input_class));
	    del=getDeBFC(T1->rightClasse(input_class),reference_class);
	    insert= getInBFC(T1->rightClasse(input_class),reference_class);
	    sub = getSuBFC(T1->rightClasse(input_class),reference_class);
	  }
	  break;
	case 2 :
	  {
	    _choicesFClass->putFirst(input_class,reference_class,T2->rightClasse(reference_class));
	    del =  getDeBFC(input_class,T2->rightClasse(reference_class));
	    
	    insert= getInBFC(input_class,T2->rightClasse(reference_class));
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
  _forestclasses[input_class+1][reference_class+1] = max;
  
}


void Matching_PO_Local::distanceBetweenClass(int input_class,int reference_class){

  int i;
  DistanceType dist =0;
  vector<int>::const_iterator begin_node;
  int ni = T1->getNbRootsInClass(input_class);
  int nj = T2->getNbRootsInClass(reference_class);
  NodeList* input_list= T1->getRootsInClass(input_class);
  NodeList* reference_list=T2->getRootsInClass(reference_class);

  LocalMatchPath _restrMapp;
  
  VertexVector _restrMappList; 
  _restrMappList.resize(ni+nj+3,EMPTY_NODE);  
  
  if(ni ==0 && nj ==0)
    dist=0;
  else{
    if (ni==0)
      {
	_restrMappList[1]=2;
	for (i=1;i<=nj;i++) { _restrMappList[i+1]=1; }
	
	dist=0;
      }
    else
      {
	if (nj==0){
	  _restrMappList[2]=1;
	  for (i=1;i<=ni;i++) { _restrMappList[i]=ni+1; }
	  dist=0;
	}
	
	else
	  {
	    _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*(_distances->getDistanceTable())); 
	    _restrMapp.make(*input_list,*reference_list);
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
	    
	    
	    dist=_restrMapp.maxCostFlow(_restrMappList);  
	    //dist = 0;
	  }
      }
  } 

  DistanceType del,insert,sub;
  del=0;insert=0;sub=0;
  
  if(ni!=0 || nj !=0){
    if (ni==0)
      {
	insert=0;
	del=0;
	sub=0;
      }
    else if (nj==0)
      {
	insert=0;
	del=0;
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
	    _choicesFClass->putLast(input_class,reference_class,_restrMapp.who(_restrMappList[i]));
	    
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


DistanceType  Matching_PO_Local::match(){
  
  int i,j;
  int size1 = T1->getNbClass();
  int size2 = T2->getNbClass();
  vector<int>::const_iterator begin_node1;
  vector<int>::const_iterator begin_node2;
  DistanceType d;  
  DistanceType D=0;
  
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
	      d=getDBT(input_root,reference_root);
		if (d>D){
		  D = d;
		  _i_v = input_root;
		  _r_v = reference_root; 
		}
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
  //_alignement->append(0,0,getSubstitutionCost(0,0));
  return D;
} 

