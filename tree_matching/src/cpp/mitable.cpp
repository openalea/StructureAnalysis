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


#include "mitable.h"

//const int EMPTY_NODE=-1;

// --------------
// Constructeur
// --------------
MatchingIntTable::MatchingIntTable(TreeGraph& input_tree, TreeGraph& reference_tree, NodeCost& node_distance)
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
  if (!T1->isNull()) inputTreeToEmpty(0);
// Si T2 est non nul on calule le cout de son insertion pour mettre a jour les vecteurs
  if (!T2->isNull()) referenceTreeFromEmpty(0);
  int r_size = T1->getDepth()*T1->getDegree();
  int c_size = T2->getNbVertex();
  int s_size = T1->getNbVertex();
  _treeDistTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestDistTable.resize(r_size, c_size, s_size);
};

void MatchingIntTable::make(TreeGraph& input_tree, TreeGraph& reference_tree, NodeCost& node_distance)
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
  //  cout<<"s_size "<<s_size<<endl;
  _treeDistTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestDistTable.resize(r_size, c_size, s_size);
// Si T1 est non nul on calule le cout de son effacement pour mettre a jour les vecteurs
  if (!T1->isNull()) inputTreeToEmpty(0);
// Si T2 est non nul on calule le cout de son insertion pour mettre a jour les vecteurs
// et on met a jour le manager d'index
  if (!T2->isNull()) 
  {
    _treeDistTable.openIntVector(T1->getNbVertex());
    referenceTreeFromEmpty(0);
  }
};

// -----------------------------------------------------------
// Calcule le cout d'effacement de l'arbre de racine vertex
//-------------------------------------------------------------
IntType MatchingIntTable::inputTreeToEmpty(int vertex)
{
// Le cout est initialise a celui de l'effacement du sommet
  IntType cost = getDCost(vertex);
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
IntType MatchingIntTable::inputForestToEmpty(int vertex)
{
  IntType cost=0;
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

// -----------------------------------------------------------
// Calcule le cout d'insertion de l'arbre de racine vertex
//-------------------------------------------------------------
IntType MatchingIntTable::referenceTreeFromEmpty(int vertex)
{
// Le cout initial est celui de l'insertion du noeud vertex
  IntType cost=getICost(vertex);
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
  _treeDistTable.putInt(cost,T1->getNbVertex(),vertex);
  return(cost);
}

// -----------------------------------------------
// Calcule le cout de l'insertion de la sous foret 
// de sommet vertex.
// -----------------------------------------------
IntType MatchingIntTable::referenceForestFromEmpty(int vertex)
{
  IntType cost = 0;
  for (int i=1;i<=T2->getNbChild(vertex);i++)
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
IntType MatchingIntTable::getDBT(int input_vertex,int reference_vertex) const 
{
  IntType distance;
// Si les deux noeud consideres sont vide alors la distance est nulle
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=0;};
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
      distance=_treeDistTable.getInt(input_vertex,reference_vertex);
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
IntType MatchingIntTable::getDBF(int input_vertex,int reference_vertex) const 
{
  IntType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) distance=0;
// Si le noeud initial est vide
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
      distance=_forestDistTable.getInt(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

// ----------------------------------------
// Insere la distance entre deux forets
// dans le tableau de distance entre forets
// -----------------------------------------
void  MatchingIntTable::putDBF(int input_vertex,int reference_vertex ,IntType new_distance)
{
  _forestDistTable.putInt(new_distance,input_vertex,reference_vertex);
}

// ----------------------------------------
// Insere la distance entre deux arbres
// dans le tableau de distance entre arbres
// -----------------------------------------
void  MatchingIntTable::putDBT(int input_vertex,int reference_vertex ,IntType new_distance)
{
  _treeDistTable.putInt(new_distance,input_vertex,reference_vertex);
}

void MatchingIntTable::openIntVector(int input_vertex)
{
  _treeDistTable.openIntVector(input_vertex);
  _treeDistTable.putInt(_inputTreeToEmpty[input_vertex],input_vertex,T2->getNbVertex());
  _forestDistTable.openIntVector(input_vertex);
}


void MatchingIntTable::closeIntVector(int input_vertex)
{
  _treeDistTable.closeIntVector(input_vertex);
  _forestDistTable.closeIntVector(input_vertex);
}

// ----------------------------------------
// Renvoie le cout d'insertion d'un noeud
// ----------------------------------------
IntType MatchingIntTable::getICost(int& vertex) const
{
  TreeNode* node=T2->getNode(vertex);
  IntType cost =1;
  
	
  return(cost);
}
// ----------------------------------------
// Renvoie le cout d'effacement d'un noeud
// ----------------------------------------
IntType MatchingIntTable::getDCost(int& vertex ) const
{
  TreeNode* node=T1->getNode(vertex);
  IntType cost=1;
    
  return(cost);
}

// ----------------------------------------------
// Renvoie le cout de substitution de deux noeuds
// ----------------------------------------------
IntType MatchingIntTable::getCCost(int vertex1 ,int vertex2) const
{
  return(0);
}



















