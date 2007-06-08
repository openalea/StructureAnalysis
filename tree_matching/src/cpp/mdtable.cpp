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


#include "mdtable.h"


// --------------
// Constructeur
// --------------
MatchingDistanceTable::MatchingDistanceTable(TreeGraph& input_tree, TreeGraph& reference_tree, NodeCost& node_distance)
{
// T1 recoit l'arbre initial
  T1 = &input_tree;
// et T2 recoit l'arbre de reference
  T2 = &reference_tree;
  ND = &node_distance;
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
};

void MatchingDistanceTable::make(TreeGraph& input_tree, TreeGraph& reference_tree, NodeCost& node_distance)
{
// T1 recoit l'arbre initial
  T1 = &input_tree;
// et T2 recoit l'arbre de reference
  T2 = &reference_tree;
  ND = &node_distance;
// le vecteur contenant les distances des noeuds de T1 au noeud vide est de taille |T1|
  _inputTreeToEmpty.resize(T1->getNbVertex());
// le vecteur contenant les distances des noeuds de T2 au noeud vide est de taille |T2|
  _referenceTreeFromEmpty.resize(T2->getNbVertex());
  int r_size = T1->getDepth()*T1->getDegree();
  int c_size = T2->getNbVertex();
  int s_size = T1->getNbVertex();
  _treeDistTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestDistTable.resize(r_size, c_size, s_size);
// Si T1 est non nul on calule le cout de son effacement pour mettre a jour les vecteurs
  if (!T1->isNull()) inputTreeToEmpty(0);
// Si T2 est non nul on calule le cout de son insertion pour mettre a jour les vecteurs
// et on met a jour le manager d'index
  if (!T2->isNull()) 
  {
    _treeDistTable.openDistanceVector(T1->getNbVertex());
    referenceTreeFromEmpty(0);
  }
};

void MatchingDistanceTable::putInputTreeToEmpty(int vertex ,DistanceType cost){
  _inputTreeToEmpty[vertex] = cost;
}
DistanceType MatchingDistanceTable::getInputTreeToEmpty(int vertex){
  return _inputTreeToEmpty[vertex];
}
void MatchingDistanceTable::putReferenceTreeFromEmpty(int vertex,DistanceType cost){
  _treeDistTable.putDistance(cost,T1->getNbVertex(),vertex);
  _referenceTreeFromEmpty[vertex]=cost;
}
DistanceType MatchingDistanceTable::getReferenceTreeFromEmpty(int vertex){
  return _referenceTreeFromEmpty[vertex];
}


// -----------------------------------------------------------
// Calcule le cout d'effacement de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType MatchingDistanceTable::inputTreeToEmpty(int vertex)
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
DistanceType MatchingDistanceTable::inputForestToEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=1;i<=T1->getNbChild(vertex);i++)
  {
// le cout d'effacement de la foret est la somme des
// cout d'effacement de tous les sous arbres dont la racine 
// est un fils du sommet
    cost=cost+inputTreeToEmpty(T1->child(vertex,i));
  }
  //  assert(cost>=0);
  return(cost);
}

// ---------------------------------------------------------
// Calcule le cout d'effacement de la foret Ordonnée issue de vertex
// ---------------------------------------------------------
DistanceType MatchingDistanceTable::inputOrderedForestToEmpty(int vertex) 
{
  DistanceType cost=MINDIST;
  int fat = T1->father(vertex);
  int ni = T1->getNbChild(fat);
  for (int i=1;i<=ni;i++)
  {
// le cout de l'insertion est la somme des couts des insertions des 
// arbres issus de vertex.
    int brother = T1->child(fat,i);
    if (brother>=vertex)
      cost=cost+inputTreeToEmpty(brother);
  }
  //  assert(cost>=0);
  return(cost);
}

// -----------------------------------------------------------
// Calcule le cout d'insertion de l'arbre de racine vertex
//-------------------------------------------------------------
DistanceType MatchingDistanceTable::referenceTreeFromEmpty(int vertex)
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
DistanceType MatchingDistanceTable::referenceForestFromEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=1;i<=T2->getNbChild(vertex);i++)
  {
// le cout de l'insertion est la somme des couts des insertions des 
// arbres issus de vertex.
    cost=cost+referenceTreeFromEmpty(T2->child(vertex,i));
  }
  //  assert(cost>=0);
  return(cost);
}


// -----------------------------------------------
// Calcule le cout de l'insertion de la sous foret 
// de sommet vertex.
// -----------------------------------------------
DistanceType MatchingDistanceTable::referenceOrderedForestFromEmpty(int vertex) 
{
  DistanceType cost=MINDIST;
  int fat = T2->father(vertex);
  int ni = T2->getNbChild(fat);
  //  int right = T2->rightBrother(vertex);
  for (int i=1;i<=ni;i++)
  {
// le cout de l'insertion est la somme des couts des insertions des 
// arbres issus de vertex.
    int brother = T2->child(fat,i);
    if (brother>=vertex)
      cost=cost+referenceTreeFromEmpty(brother);
  }
  //  assert(cost>=0);
  return(cost);
}
// ------------------------------------------------
// Renvoie la distance entre deux arbres de racines
// input_vertex et reference_vertex
// ------------------------------------------------
DistanceType MatchingDistanceTable::getDBT(int input_vertex,int reference_vertex)  const
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
DistanceType MatchingDistanceTable::getDBF(int input_vertex,int reference_vertex)  const
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

// ------------------------------------------------------
// Renvoie la distance entre deux forets ordonnées dont 
// le premier arbre a pour racines
// input_vertex et reference_vertex le raisonnement est 
// a celui qui est vu dans la fonction precedente
// ------------------------------------------------------
DistanceType MatchingDistanceTable::getDBOrderedF(int input_vertex,int reference_vertex) 
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;}
// Si le noeud initial est vide
  else{
    if (input_vertex==EMPTY_NODE)
      {
	distance=referenceOrderedForestFromEmpty(reference_vertex);
      }
    else
      {
	if (reference_vertex==EMPTY_NODE)
	  {
	    distance=inputOrderedForestToEmpty(input_vertex);
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
void  MatchingDistanceTable::putDBF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _forestDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

// ----------------------------------------
// Insere la distance entre deux arbres
// dans le tableau de distance entre arbres
// -----------------------------------------
void  MatchingDistanceTable::putDBT(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
  _treeDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void MatchingDistanceTable::openDistancesVector(int input_vertex)
{
  _treeDistTable.openDistanceVector(input_vertex);
  _treeDistTable.putDistance(_inputTreeToEmpty[input_vertex],input_vertex,T2->getNbVertex());
  _forestDistTable.openDistanceVector(input_vertex);
}


void MatchingDistanceTable::closeDistancesVector(int input_vertex)
{
  _treeDistTable.closeDistanceVector(input_vertex);
  _forestDistTable.closeDistanceVector(input_vertex);
}

// ----------------------------------------
// Renvoie le cout d'insertion d'un noeud
// ----------------------------------------
DistanceType MatchingDistanceTable::getICost(int& vertex) const
{
  TreeNode* node=T2->getNode(vertex);
  DistanceType cost;
  
  switch(ND->type())
  { 
    case TOPOLOGIC  : cost=ND->getInsertionCost(node);break;
    case LOCAL_TOPO : cost=((TopologicalLocalNodeCost*) ND)->getInsertionCost(node);break;
    case WEIGTH     : cost=((WeightedNodeCost*) ND)->getInsertionCost(node);break;
	case SCORE     : cost=((WeightedNodeCost*) ND)->getInsertionCost(node);break;
    case MATRIX     : cost=((MatrixNodeCost*)   ND)->getInsertionCost(node);break;
    default         : assert(0);break;
  }
	
  return(cost);
}
// ----------------------------------------
// Renvoie le cout d'effacement d'un noeud
// ----------------------------------------
DistanceType MatchingDistanceTable::getDCost(int& vertex ) const
{
  TreeNode* node=T1->getNode(vertex);
  DistanceType cost;
  
  switch(ND->type())
  { 
     case TOPOLOGIC  : { cost = ND->getDeletionCost(node); };break;
     case LOCAL_TOPO : { cost = ((TopologicalLocalNodeCost*) ND)->getDeletionCost(node); };break;
     case WEIGTH     : { cost = ((WeightedNodeCost*) ND)->getDeletionCost(node); };break;
	 case SCORE     : cost=((WeightedNodeCost*) ND)->getDeletionCost(node);break;
     case MATRIX     : { cost = ((MatrixNodeCost*) ND)->getDeletionCost(node); };break;
     default         : { assert(0); };break;
  }
  
  return(cost);
}

// ----------------------------------------------
// Renvoie le cout de substitution de deux noeuds
// ----------------------------------------------
DistanceType MatchingDistanceTable::getCCost(int vertex1 ,int vertex2) const
{
// node1 et node2 recoivent les noeuds correspondant aux indices
	TreeNode* node1=T1->getNode(vertex1);
	TreeNode* node2=T2->getNode(vertex2);
	DistanceType cost;

	switch(ND->type())
       	{ 
// le cout est donne par la classe nodecost et le node cost ND
	case TOPOLOGIC  : {cost=ND->getChangingCost(node1,node2);};break;
	case LOCAL_TOPO : {cost=((TopologicalLocalNodeCost*) ND)->getChangingCost(node1,node2);};break;
	case WEIGTH     : {cost=((WeightedNodeCost*) ND)->getChangingCost(node1,node2);};break;
	case SCORE     : cost=((WeightedNodeCost*) ND)->getChangingCost(node1,node2);break;
	case MATRIX     : {cost=((MatrixNodeCost*) ND)->getChangingCost(node1,node2);};break;
	default         : assert(0);break;
        }
	
	return(cost);
}
