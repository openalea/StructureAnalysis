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


#include "newmatching.h"


NewMatching::NewMatching(TreeGraph& imput,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&imput;
  T2=&reference;
  ND=&nodeDistance;
  _distances.make(*T1,*T2,nodeDistance);
  // _choices est un tableau de listes retenant les tentatives successives alignements durant l'algo.
  // c'est donc un tableau de |T1| lignes et |T2| colonnes initialise a 0
  _choices.resize(T1->getNbVertex(),T2->getNbVertex());
  // constante qui va permettre de calculer l'alignement restreint
  _restrMapp.link(I_MAX(T1->getDegree(),T2->getDegree()),*_distances.getDistanceTable());
}
// -------------
// Destructeur
// -------------
NewMatching::~NewMatching()
{
}

 


// ----------------------------------------------------------------------------------
// Calcule la distance entre les deux arbres T1[imput_vertex] et T2[reference_vertex]
// ----------------------------------------------------------------------------------
DistanceType NewMatching::distanceBetweenTree(int imput_vertex,int reference_vertex)
{
  // On stocke dans ni et nj le nombre d'enfants de T1[i] et T2[j]
  int ni=T1->getNbChild(imput_vertex);
  int nj=T2->getNbChild(reference_vertex);

  DistanceType cost3;
  DistanceType MIN=2*MAXDIST;

  //----------------------------------------------------------------------------------
  //Case3 : We evaluate the matching between the imput_forest and the reference_forest
  // On evalue la mise en correspondance des arbres des deux forets issues de T1 et T2
  //----------------------------------------------------------------------------------
  // Le cout est celui de l'alignement des deux forets
  // plus celui de l'echange de T1(i) en T2(j)
  cost3=getDBF(imput_vertex,reference_vertex);
  cost3=cost3+_distances.getCCost(imput_vertex,reference_vertex);
  // On conserve le cout s'il est inferieur au precedent
  if (cost3<MIN)  MIN=cost3;

  //-----------------------------------
  // We maintain the matching lists
  // mise a jour des listes d'alignement
  //-----------------------------------
  _choices.putFirst(imput_vertex,reference_vertex,-1);
  _choices.putLast(imput_vertex,reference_vertex,reference_vertex);
  _choices.putFirst(imput_vertex,reference_vertex,3);
  // On range dans le tableau des distances, la distance entre les arbres de racines imput_vertex et
  // reference_vertex.
  _distances.putDBT(imput_vertex,reference_vertex,MIN);
  return(MIN);

}



// -------------------------------
// Calcule la distance entre deux
// forets
// -------------------------------
DistanceType NewMatching::distanceBetweenForest(int imput_vertex,int reference_vertex)
{
  // ni et nj representent le nombre de forets a comparees
  int ni=T1->getNbChild(imput_vertex);
  int nj=T2->getNbChild(reference_vertex);
  DistanceType cost3;
  DistanceType DIST;
  DIST=MAXDIST;


  //---------------------------------------------------------------------------------------------
  //Case 3 : We evaluate the restricted mapping between the imput_forest and the reference_forest
  // On evalue l'alignement restreint entre les deux forets
  //---------------------------------------------------------------------------------------------

  // On fabrique le graphe de flot necessaire a la resolution du probleme
  NodeList* input_list=new NodeList();
  NodeList* reference_list=new NodeList();
  for (int s1=1;s1<=ni;s1++) { input_list->push_back(T1->child(imput_vertex,s1)); };
  for (int s2=1;s2<=nj;s2++) { reference_list->push_back(T2->child(reference_vertex,s2)); };
  _restrMapp.make(*input_list,*reference_list);
  _restrMappList.resize(ni+nj+3,EMPTY_NODE);


  // THE IMPUT FOREST IS EMPTY_TREE
  // All the reference vertices are paired with empty
  // Si la foret initiale est vide, il faut inserer toutes les arbres de
  // la foret de reference et tous les noeuds de references sont associes 
  // avec le noeud vide
  if (ni==0) 
    { 
      _restrMappList[1]=2;
      for (int  i=1;i<=nj;i++) { _restrMappList[i+1]=1; };
      cost3=getDBF(EMPTY_TREE,reference_vertex);
    }
  else
    {
      // THE REFERENCE FOREST IS EMPTY_TREE
      // All the imput vertices are paired with empty
      // Si c'est l'arbre de reference qui est vide,
      // il faut supprimer la foret initiale et tous les 
      // noeuds de cette foret seront associer avec un
      // noeud vide
      if (nj==0) 
	{ 
	  _restrMappList[2]=1;
	  for (int i=1;i<=ni;i++) { _restrMappList[i]=ni+1; };
	  cost3=getDBF(imput_vertex,EMPTY_TREE); 
	}
      else
	{
	  //BOTH FOREST ARE NOT EMPTY_TREE
	  // A retricted mapping must be calculated
	  // Sinon on resout le probleme de flot maximum de cout minimum		  
	  cost3=_restrMapp.minCostFlow(_restrMappList);
		 
	}
    }

	
  if (cost3<=DIST)  DIST=cost3;

  //-------------------------------
  //We maintain the matching lists
  // On maintient les listes d'alignement
  //-------------------------------
  _choices.createList(imput_vertex,reference_vertex);
  _choices.putFirst(imput_vertex,reference_vertex,3);


  for (int i=1;i<=T1->getNbChild(imput_vertex);i++)
    {
      _choices.putLast(imput_vertex,reference_vertex,_restrMapp.who(_restrMappList[i]));
    }
	
  delete (NodeList*) input_list;
  delete (NodeList*) reference_list;
  _distances.putDBF(imput_vertex,reference_vertex,DIST);
  return(DIST);

}


