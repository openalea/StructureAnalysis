/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source$
 *       $Id: mdtable.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include "mdtable.h"


// --------------
// Constructeur
// --------------
MatchingDistanceTable::MatchingDistanceTable(const TreeGraphPtr& input_tree, 
											 const TreeGraphPtr& reference_tree, 
											 const NodeCostPtr node_distance):
  T1(input_tree),T2(reference_tree), ND(node_distance)
{
  assert(T1);
  assert(T2);
  assert(ND);

  _n1 = T1->getNbVertex();
  _n2 = T2->getNbVertex();

  _inputTreeToEmpty.resize(_n1);
  _referenceTreeFromEmpty.resize(_n2);
  _treeDistTable.resize(_n1+1);
  _forestDistTable.resize(_n1+1);
  for (int i = 0; i<=_n1 ; i++){
    _treeDistTable[i].resize(_n2+1);
    _forestDistTable[i].resize(_n2+1);
  }
}

void MatchingDistanceTable::putInputTreeToEmpty(int vertex ,DistanceType cost){
  //_inputTreeToEmpty[vertex] = cost;
  _treeDistTable[vertex][_n2] = cost;
}

DistanceType MatchingDistanceTable::getInputTreeToEmpty(int vertex){
  //  return _inputTreeToEmpty[vertex];
  return _treeDistTable[vertex][_n2];
}
void MatchingDistanceTable::putReferenceTreeFromEmpty(int vertex,DistanceType cost){
  //_referenceTreeFromEmpty[vertex]=cost;
  _treeDistTable[_n1][vertex] = cost;
}
DistanceType MatchingDistanceTable::getReferenceTreeFromEmpty(int vertex){
  //  return _referenceTreeFromEmpty[vertex];
  return _treeDistTable[_n1][vertex];
}


// -----------------------------------------------------------
// Calcule le cout d'effacement de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType MatchingDistanceTable::inputTreeToEmpty(int vertex)
{
//   DistanceType cost = getDCost(vertex);
//   if (!T1->isLeaf(vertex))
//     cost+=inputForestToEmpty(vertex);
//   _inputTreeToEmpty[vertex] = cost;
//   return(cost);
  DistanceType cost = getDCost(vertex);
  if (!T1->isLeaf(vertex))
    cost+=_forestDistTable[vertex][_n2];
  _treeDistTable[vertex][_n2] = cost;
  return(cost);
}

// ---------------------------------------------------------
// Calcule le cout d'effacement de la foret issue de vertex
// le cout d'effacement de la foret est la somme des
// cout d'effacement de tous les sous arbres dont la racine 
// est un fils du sommet
// ---------------------------------------------------------
DistanceType MatchingDistanceTable::inputForestToEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=0;i<T1->getNbChild(vertex);i++){

//     cost+=inputTreeToEmpty(T1->child(vertex,i));
    cost+=_treeDistTable[T1->child(vertex,i)][_n2];
  }
  _forestDistTable[vertex][_n2] = cost;
  return(cost);
}


// -----------------------------------------------------------
// Calcule le cout d'insertion de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType MatchingDistanceTable::referenceTreeFromEmpty(int vertex)
{
  DistanceType cost = getICost(vertex);
  if (!T2->isLeaf(vertex))
    cost+=_forestDistTable[_n1][vertex];
  _treeDistTable[_n1][vertex] = cost;
  return(cost);

}

// -----------------------------------------------
// Calcule le cout de l'insertion de la sous foret 
// de sommet vertex.
// -----------------------------------------------
DistanceType MatchingDistanceTable::referenceForestFromEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=0;i<T2->getNbChild(vertex);i++)
    cost+=_treeDistTable[_n1][T2->child(vertex,i)];
  _forestDistTable[_n1][vertex] = cost;
  return(cost);
}



// ------------------------------------------------
// Renvoie la distance entre deux arbres de racines
// input_vertex et reference_vertex
// ------------------------------------------------
DistanceType MatchingDistanceTable::getDBT(int input_vertex,int reference_vertex)  const
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) 
    distance = MINDIST;
  if (input_vertex==EMPTY_NODE)
//     distance=_referenceTreeFromEmpty[reference_vertex];
    distance = _treeDistTable[_n1][reference_vertex];
  else{
    if (reference_vertex==EMPTY_NODE)
//      distance=_inputTreeToEmpty[input_vertex];
      distance = _treeDistTable[input_vertex][_n2];
    else
      distance=_treeDistTable[input_vertex][reference_vertex];
  }
  return(distance);
}

// ------------------------------------------------------
// Renvoie la distance entre deux forets de racines
// input_vertex et reference_vertex le raisonnement est 
// a celui qui est vu dans la fonction precedente
// ------------------------------------------------------
DistanceType MatchingDistanceTable::getDBF(int input_vertex,int reference_vertex)  const
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) 
    distance=MINDIST;
  else{
    if (input_vertex==EMPTY_NODE)
// 	distance=_referenceTreeFromEmpty[reference_vertex]-getICost(ref_n2+1erence_vertex);
      distance = _forestDistTable[_n1][reference_vertex];
    else{
	if (reference_vertex==EMPTY_NODE)
// 	    distance=_inputTreeToEmpty[input_vertex]-getDCost(input_vertex);
	  distance = _forestDistTable[input_vertex][_n2];
	else
	    distance=_forestDistTable[input_vertex][reference_vertex];
    }
  }
  return(distance);
}


// ----------------------------------------
// Insere la distance entre deux forets
// dans le tableau de distance entre forets
// -----------------------------------------
void  MatchingDistanceTable::putDBF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _forestDistTable[input_vertex][reference_vertex]=new_distance;
}

// ----------------------------------------
// Insere la distance entre deux arbres
// dans le tableau de distance entre arbres
// -----------------------------------------
void  MatchingDistanceTable::putDBT(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _treeDistTable[input_vertex][reference_vertex] = new_distance;
}

// ----------------------------------------
// Renvoie le cout d'insertion d'un noeud
// ----------------------------------------
DistanceType MatchingDistanceTable::getICost(int& vertex) const
{  
	return ND->getInsertionCost(T2->getNode(vertex)); 
}

// ----------------------------------------
// Renvoie le cout d'effacement d'un noeud
// ----------------------------------------
DistanceType MatchingDistanceTable::getDCost(int& vertex ) const
{  
	return ND->getDeletionCost(T1->getNode(vertex)); 
}

// ----------------------------------------------
// Renvoie le cout de substitution de deux noeuds
// ----------------------------------------------
DistanceType MatchingDistanceTable::getCCost(int vertex1 ,int vertex2) const
{  return ND->getChangingCost(T1->getNode(vertex1),T2->getNode(vertex2)); }
