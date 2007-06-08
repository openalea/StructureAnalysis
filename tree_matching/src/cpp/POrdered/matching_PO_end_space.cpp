/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): A.ouangraoua (ouangrao@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/POrdered/matching_PO_end_space.cpp,v $
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


#include "matching_PO_end_space.h"


const int EMPTY_CLASS=-1;


// -------------
  // Constructeur
  // -------------
Matching_PO_End_Space::Matching_PO_End_Space(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance):
  Matching_PO(input,reference,nodeDistance)
{
}



DistanceType  Matching_PO_End_Space::match(){

  int i,j;
  int size1 = T1->getNbClass();
  int size2 = T2->getNbClass();
  vector<int>::const_iterator begin_node1;
  vector<int>::const_iterator begin_node2;
  DistanceType d;  
  DistanceType D=MAXDIST;
  
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
	      if ((input_root == 0) || (reference_root == 0)){
		if (d<D){
		  D = d;
		  _i_v = input_root;
		  _r_v = reference_root; 
		}
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
  _alignement->append(0,0,getSubstitutionCost(0,0));
  return D;
} 

