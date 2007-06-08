/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): P.ferraro (pascal.ferraro@cirad.fr)
 *
 *       $Source$
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


#include "ferraro_hanna.h"

FerraroHanna::FerraroHanna(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _distances.make(*T1,*T2,nodeDistance);
//  _insertCost.make(*T1,*T2,nodeDistance);
//  _deleteCost.make(*T1,*T2,nodeDistance);
//  _substitutionCost.make(*T1,*T2,nodeDistance);
//  _sumNbCaseVector = CaseVector(3,0);
//  _nbCaseVector = CaseVector(7,0);
  _distanceMatrix = DistanceVectorTable(T1->getNbVertex());
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}


DistanceType FerraroHanna::distanceBetweenTree(int input_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(input_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType MIN;
//  DistanceType insert,del,sub;
  
  int child_input,child_ref,last_child_input,last_child_ref;
  child_input = T1->child(input_vertex,1);
  child_ref = T2->child(reference_vertex,1);
  last_child_input = T1->child(input_vertex,ni);
  last_child_ref = T2->child(reference_vertex,nj);
  
  if (((T1->getNode(input_vertex))->depth()!=2)||((T2->getNode(reference_vertex))->depth()!=2)){
	MIN = getDBF(child_input,child_ref)+_distances.getCCost(input_vertex,reference_vertex);
  }
  else{
//	  cerr<<T1->father(input_vertex)<<" - "<<T2->father(reference_vertex)<<endl;
	MIN = getPermutation(input_vertex,reference_vertex)+_distances.getCCost(input_vertex,reference_vertex);
  }
//  cerr<<" c_cost = "<<_distances.getCCost(input_vertex,reference_vertex)<<endl;

  _choices.createList(input_vertex,reference_vertex);
  _choices.putFirst(input_vertex,reference_vertex,3);
  _choices.putLast(input_vertex,reference_vertex,reference_vertex);
//  insert=getInBF(child_input,child_ref);
//  del=getDeBF(child_input,child_ref);
//  sub=getSuBF(child_input,child_ref)+_distances.getCCost(input_vertex,reference_vertex);

  // On range dans le tableau des distances, la distance entre les arbres de racines input_vertex et
  // reference_vertex.
  //  _choices.putFirst(input_vertex,reference_vertex,MTC);
//  _insertCost.putDBT(input_vertex,reference_vertex,insert);
//  _deleteCost.putDBT(input_vertex,reference_vertex,del);
//  _substitutionCost.putDBT(input_vertex,reference_vertex,sub);
  
  _distances.putDBT(input_vertex,reference_vertex,MIN);
  return(MIN);
}

// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType FerraroHanna::distanceBetweenForest(int input_vertex,int reference_vertex)
{
// ni et nj representent le nombre de forets a comparees
  DistanceType cost1,cost2,cost3;
//DistanceType insert,del,sub;
  DistanceType DIST;
  int MFC=0;

  DIST=MAXDIST;

  	int input_right_brother=T1->rightBrother(input_vertex);

	int reference_right_brother = T2->rightBrother(reference_vertex);
 

	cost1=getDBF(input_right_brother,reference_vertex)+getDBT(input_vertex,EMPTY_TREE);
	if (cost1<=DIST) {
		DIST=cost1;
		MFC=1;
	}
 
	cost2=getDBF(input_vertex,reference_right_brother)+getDBT(EMPTY_TREE,reference_vertex);

	if (cost2<=DIST) {	
		DIST=cost2;
		MFC=2;
	}
	cost3=getDBF(input_right_brother,reference_right_brother)+getDBT(input_vertex,reference_vertex);

	if (cost3<=DIST) { 
		DIST=cost3; 
		MFC=3; 
	}

	//-------------------------------
	//We maintain the matching lists
	// On maintient les listes d'alignement
	//-------------------------------
	//   _choices.createList(input_vertex,reference_vertex);
	_choices.putLast(input_vertex,reference_vertex,MFC);
  
  
  switch(MFC)
    {
    case 1 :
      {
	_choices.putLast(input_vertex,reference_vertex,input_right_brother);
      }
      break;
    case 2 :
      {
	_choices.putLast(input_vertex,reference_vertex,reference_right_brother);
      }
      break;
    case 3 :{
			_choices.putLast(input_vertex,reference_vertex,reference_vertex);
      }
      break;
    default :   break;
    }
  
  _distances.putDBF(input_vertex,reference_vertex,DIST);
  return(DIST);
}

DistanceType FerraroHanna::getPermutation(int input_vertex,int reference_vertex){
	  //---------------------------------------------------------------------------------------------
     //Case 3 : We evaluate the restricted mapping between the input_forest and the reference_forest
     // On evalue l'alignement restreint entre les deux forets
     //---------------------------------------------------------------------------------------------
     // On fabrique le graphe de flot necessaire a la resolution du probleme
	DistanceType DIST;
    NodeList input_list= NodeList();
    NodeList reference_list= NodeList();
    int ni=T1->getNbChild(input_vertex);
    int nj=T2->getNbChild(reference_vertex);

    for (int s1=1;s1<=ni;s1++) { input_list.push_back(T1->child(input_vertex,s1)); }
    for (int s2=1;s2<=nj;s2++) { reference_list.push_back(T2->child(reference_vertex,s2)); }
    _restrMapp.make(input_list,reference_list);
    _restrMappList.resize(ni+nj+3,EMPTY_NODE);
    // THE INPUT FOREST IS EMPTY_TREE
    // All the reference vertices are paired with empty
    // Si la foret initiale est vide, il faut inserer toutes les arbres de
    // la foret de reference et tous les noeuds de references sont associes
    // avec le noeud vide

          //BOTH FOREST ARE NOT EMPTY_TREE
          // A retricted mapping must be calculated
          // Sinon on resout le probleme de flot maximum de cout minimum
          DIST=_restrMapp.minCostFlow(_restrMappList);
    return DIST;
  }
