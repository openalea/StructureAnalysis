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
StdMatchingDistanceTable::StdMatchingDistanceTable(const TreeGraphPtr& input_tree, 
						   const TreeGraphPtr& reference_tree, 
						   const NodeCostPtr node_distance):
  T1(input_tree),T2(reference_tree), ND(node_distance)
{
  _type = STD;
  assert(T1);
  assert(T2);
  assert(ND);

  _n1 = T1->getNbVertex();
  _n2 = T2->getNbVertex();

  //_inputTreeToEmpty.resize(_n1);
  //_referenceTreeFromEmpty.resize(_n2);
  _treeDistTable.resize(_n1+1);
  _forestDistTable.resize(_n1+1);
  for (int i = 0; i<=_n1 ; i++){
    _treeDistTable[i].resize(_n2+1);
    _forestDistTable[i].resize(_n2+1);
  }
}

void StdMatchingDistanceTable::putInputTreeToEmpty(int vertex ,DistanceType cost){
  //_inputTreeToEmpty[vertex] = cost;
  _treeDistTable[vertex][_n2] = cost;
}

DistanceType StdMatchingDistanceTable::getInputTreeToEmpty(int vertex){
  //  return _inputTreeToEmpty[vertex];
  return _treeDistTable[vertex][_n2];
}
void StdMatchingDistanceTable::putReferenceTreeFromEmpty(int vertex,DistanceType cost){
  //_referenceTreeFromEmpty[vertex]=cost;
  _treeDistTable[_n1][vertex] = cost;
}
DistanceType StdMatchingDistanceTable::getReferenceTreeFromEmpty(int vertex){
  //  return _referenceTreeFromEmpty[vertex];
  return _treeDistTable[_n1][vertex];
}


// -----------------------------------------------------------
// Calcule le cout d'effacement de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType StdMatchingDistanceTable::inputTreeToEmpty(int vertex)
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
DistanceType StdMatchingDistanceTable::inputForestToEmpty(int vertex)
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
DistanceType StdMatchingDistanceTable::referenceTreeFromEmpty(int vertex)
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
DistanceType StdMatchingDistanceTable::referenceForestFromEmpty(int vertex)
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
DistanceType StdMatchingDistanceTable::getDBT(int input_vertex,int reference_vertex)  const
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
DistanceType StdMatchingDistanceTable::getDBF(int input_vertex,int reference_vertex)  const
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
void  StdMatchingDistanceTable::putDBF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _forestDistTable[input_vertex][reference_vertex]=new_distance;
}

// ----------------------------------------
// Insere la distance entre deux arbres
// dans le tableau de distance entre arbres
// -----------------------------------------
void  StdMatchingDistanceTable::putDBT(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _treeDistTable[input_vertex][reference_vertex] = new_distance;
}

// ----------------------------------------
// Renvoie le cout d'insertion d'un noeud
// ----------------------------------------
DistanceType StdMatchingDistanceTable::getICost(int& vertex) const
{  
 return ND->getInsertionCost(T2->getNode(vertex)); 
}

// ----------------------------------------
// Renvoie le cout d'effacement d'un noeud
// ----------------------------------------
DistanceType StdMatchingDistanceTable::getDCost(int& vertex ) const
{  
	return ND->getDeletionCost(T1->getNode(vertex)); 
}

// ----------------------------------------------
// Renvoie le cout de substitution de deux noeuds
// ----------------------------------------------
DistanceType StdMatchingDistanceTable::getCCost(int vertex1 ,int vertex2) const
{  
  return ND->getChangingCost(T1->getNode(vertex1),T2->getNode(vertex2)); 
}

// ----------------------------------------------
// Renvoie le cout de Merging de deux noeuds
// ----------------------------------------------
DistanceType StdMatchingDistanceTable::getMCost(vector<int> vertex1 ,int vertex2) const
{  
  vector<TreeNodePtr> path;
  for (int i = 0; i<vertex1.size();i++){
    path.push_back( T1->getNode(vertex1[i]));
  }
  return ND->getMergingCost(path,T2->getNode(vertex2)); 
}

DistanceType StdMatchingDistanceTable::getSCost(int vertex1 ,vector<int> vertex2) const
{  
  vector<TreeNodePtr> path;
  for (int i = 0; i<vertex2.size();i++)
    path.push_back( T2->getNode(vertex2[i]));
  return ND->getSplittingCost(T1->getNode(vertex1),path); 
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


// --------------
// Constructeur
// --------------
CompactMatchingDistanceTable::CompactMatchingDistanceTable(const TreeGraphPtr& input_tree, 
						   const TreeGraphPtr& reference_tree, 
						   const NodeCostPtr node_distance):
  T1(input_tree),T2(reference_tree), ND(node_distance)
{
  _type = COMPACT;
  assert(T1);
  assert(T2);
  assert(ND);

  _n1 = T1->getNbVertex();
  _n2 = T2->getNbVertex();
// le vecteur contenant les distances des noeuds de T1 au noeud vide est de taille |T1|
  _inputTreeToEmpty.resize(T1->getNbVertex());
// le vecteur contenant les distances des noeuds de T2 au noeud vide est de taille |T2|
  _referenceTreeFromEmpty.resize(T2->getNbVertex());
// Si T1 est non nul on calule le cout de son effacement pour mettre a jour les vecteurs
  //if (!T1->isNull()) inputTreeToEmpty(0);
// Si T2 est non nul on calule le cout de son insertion pour mettre a jour les vecteurs
  //if (!T2->isNull()) referenceTreeFromEmpty(0);
  int r_size = T1->getDepth()*T1->getDegree();
  int c_size = T2->getNbVertex();
  int s_size = T1->getNbVertex();
  r_size = s_size;
  _treeDistTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestDistTable.resize(1+r_size, 1+c_size, 1+s_size);
  if (!T1->isNull()) inputTreeToEmpty(0);
// Si T2 est non nul on calule le cout de son insertion pour mettre a jour les vecteurs
// et on met a jour le manager d'index
  if (!T2->isNull()) 
  {
    _treeDistTable.openDistanceVector(_n1);
    referenceTreeFromEmpty(0);
  }
}


void CompactMatchingDistanceTable::putInputTreeToEmpty(int vertex ,DistanceType cost){
  _inputTreeToEmpty[vertex] = cost;
}
DistanceType CompactMatchingDistanceTable::getInputTreeToEmpty(int vertex){
  return _inputTreeToEmpty[vertex];
}
void CompactMatchingDistanceTable::putReferenceTreeFromEmpty(int vertex,DistanceType cost){
  _treeDistTable.putDistance(cost,T1->getNbVertex(),vertex);
  _referenceTreeFromEmpty[vertex]=cost;
}
DistanceType CompactMatchingDistanceTable::getReferenceTreeFromEmpty(int vertex){
  return _referenceTreeFromEmpty[vertex];
}


// -----------------------------------------------------------
// Calcule le cout d'effacement de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType CompactMatchingDistanceTable::inputTreeToEmpty(int vertex)
{
// Le cout est initialise a celui de l'effacement du sommet
  DistanceType cost = getDCost(vertex);
// Si vertex n'est pas une feuille
  if (!T1->isLeaf(vertex))
  {
// alors le cout d'effacement de cet arbre est celui de l'effacement de
// sa racine plus celui de la foret compose des sous arbres du vertex
// racine.
    cost=cost+inputForestToEmpty(vertex);
  }
  //  assert(cost>=0);
// On place la distance calculee dans le vecteur approprie
  _inputTreeToEmpty[vertex] = cost;
  return(cost);
}

// ---------------------------------------------------------
// Calcule le cout d'effacement de la foret issue de vertex
// ---------------------------------------------------------
DistanceType CompactMatchingDistanceTable::inputForestToEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=0;i<T1->getNbChild(vertex);i++)
  {
 // le cout d'effacement de la foret est la somme des
// cout d'effacement de tous les sous arbres dont la racine 
// est un fils du sommet
    cost=cost+inputTreeToEmpty(T1->child(vertex,i));
  }
  //  assert(cost>=0);
  return(cost);
}


// -----------------------------------------------------------
// Calcule le cout d'insertion de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType CompactMatchingDistanceTable::referenceTreeFromEmpty(int vertex)
{
// Le cout initial est celui de l'insertion du noeud vertex
  DistanceType cost=getICost(vertex);
// Si l'arbre a insere n'est pas reduit a une feuille
  if (!T2->isLeaf(vertex))
  {
// Le cout est celui de l'insertion du noeud plus celui de l'insertion
// de tous les sous arbres de racine vertex.
    cost=cost+referenceForestFromEmpty(vertex);
  }
  //  assert(cost>=0);
// On insere ce cout dans le vecteur correspondant.
  _referenceTreeFromEmpty[vertex]=cost;
//
  _treeDistTable.putDistance(cost,T1->getNbVertex(),vertex);
  return(cost);
}

// -----------------------------------------------
// Calcule le cout de l'insertion de la sous foret 
// de sommet vertex.
// -----------------------------------------------
DistanceType CompactMatchingDistanceTable::referenceForestFromEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=0;i<T2->getNbChild(vertex);i++)
  {
// le cout de l'insertion est la somme des couts des insertions des 
// arbres issus de vertex.
    cost=cost+referenceTreeFromEmpty(T2->child(vertex,i));
  }
  //  assert(cost>=0);
  return(cost);
}


// ------------------------------------------------
// Renvoie la distance entre deux arbres de racines
// input_vertex et reference_vertex
// ------------------------------------------------
DistanceType CompactMatchingDistanceTable::getDBT(int input_vertex,int reference_vertex)  const
{
  DistanceType distance;
// Si les deux noeud consideres sont vide alors la distance est nulle
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
// Si le sommet initial est vide alors on a deja caluler le cout de l'insertion
// de l'arbre de sommet reference_vertex qui est stocke dans le vecteur 
// _referenceTreeFromEmpty
  if (input_vertex==EMPTY_NODE)
  {
    distance=_referenceTreeFromEmpty[reference_vertex];
  }
  else
  {
// Sinon si c'est le noeud e reference qui est vide, on egalement stocke le calcul
// de l'effacement dans _inputTreeToEmpty
    if (reference_vertex==EMPTY_NODE)
    {
      distance=_inputTreeToEmpty[input_vertex];
    }
    else
    {
// Sinon, enfin, on renvoie la distance entre ces deux arbres contenue dans le
// tableau des distances.
      distance=_treeDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
// on renvoie la distance trouvee
  return(distance);
}

// ------------------------------------------------------
// Renvoie la distance entre deux forets de racines
// input_vertex et reference_vertex le raisonnement est 
// a celui qui est vu dans la fonction precedente
// ------------------------------------------------------
DistanceType CompactMatchingDistanceTable::getDBF(int input_vertex,int reference_vertex)  const
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;}
// Si le noeud initial est vide
  else{
    if (input_vertex==EMPTY_NODE)
      {
	// alors la distance est celle de l'insertion de l'arbre de racine reference_vertex
	// moins l'insertion du noeud puisqu'on ne calcule la distance entre forets issues
	//de deux neouds
	distance=_referenceTreeFromEmpty[reference_vertex]-getICost(reference_vertex);
      }
    else
      {
	if (reference_vertex==EMPTY_NODE)
	  {
	    // si le noeude de reference est nul alors c'est le cas symetrique
	    distance=_inputTreeToEmpty[input_vertex]-getDCost(input_vertex);
	  }
	else
	  // Sinon la distance est stocke dans le tableau des distances.
	  {
	    distance=_forestDistTable.getDistance(input_vertex,reference_vertex);
	  }
      }
  }
  return(distance);
}


// ----------------------------------------
// Insere la distance entre deux forets
// dans le tableau de distance entre forets
// -----------------------------------------
void  CompactMatchingDistanceTable::putDBF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _forestDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

// ----------------------------------------
// Insere la distance entre deux arbres
// dans le tableau de distance entre arbres
// -----------------------------------------
void  CompactMatchingDistanceTable::putDBT(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _treeDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void CompactMatchingDistanceTable::openDistancesVector(int input_vertex)
{
  _treeDistTable.openDistanceVector(input_vertex);
  _treeDistTable.putDistance(_inputTreeToEmpty[input_vertex],input_vertex,T2->getNbVertex());
  _forestDistTable.openDistanceVector(input_vertex);
}


void CompactMatchingDistanceTable::closeDistancesVector(int input_vertex)
{
  _treeDistTable.closeDistanceVector(input_vertex);
  _forestDistTable.closeDistanceVector(input_vertex);
}

// ----------------------------------------
// Renvoie le cout d'insertion d'un noeud
// ----------------------------------------
DistanceType CompactMatchingDistanceTable::getICost(int& vertex) const
{
 return ND->getInsertionCost(T2->getNode(vertex)); 

}
// ----------------------------------------
// Renvoie le cout d'effacement d'un noeud
// ----------------------------------------
DistanceType CompactMatchingDistanceTable::getDCost(int& vertex ) const
{
  return ND->getDeletionCost(T1->getNode(vertex)); 
}


// ----------------------------------------------
// Renvoie le cout de substitution de deux noeuds
// ----------------------------------------------
DistanceType CompactMatchingDistanceTable::getCCost(int vertex1 ,int vertex2) const
{
  return ND->getChangingCost(T1->getNode(vertex1),T2->getNode(vertex2)); 
}

// ----------------------------------------------
// Renvoie le cout de Merging de deux noeuds
// ----------------------------------------------
DistanceType CompactMatchingDistanceTable::getMCost(vector<int> vertex1 ,int vertex2) const
{  
  vector<TreeNodePtr> path;
  for (int i = 0; i<vertex1.size();i++)
    path.push_back( T1->getNode(vertex1[i]));
  return ND->getMergingCost(path,T2->getNode(vertex2)); 
}

DistanceType CompactMatchingDistanceTable::getSCost(int vertex1 ,vector<int> vertex2) const
{  
  vector<TreeNodePtr> path;
  for (int i = 0; i<vertex2.size();i++)
    path.push_back( T2->getNode(vertex2[i]));
  return ND->getSplittingCost(T1->getNode(vertex1),path); 
}

