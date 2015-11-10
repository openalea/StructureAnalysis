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


#include "mctable.h"


// --------------
// Constructeur
// --------------
MatchingConnectedTable::MatchingConnectedTable(TreeGraph& input_tree, TreeGraph& reference_tree)
{
// T1 recoit l'arbre initial
  T1 = &input_tree;
// et T2 recoit l'arbre de reference
  T2 = &reference_tree;
  int r_size = T1->getDepth()*T1->getDegree();
  int c_size = T2->getNbVertex();
  int s_size = T1->getNbVertex();
  _treeConnectedTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestConnectedTable.resize(r_size, c_size, s_size);
}

void MatchingConnectedTable::make(TreeGraph& input_tree, TreeGraph& reference_tree)
{
// T1 recoit l'arbre initial
  T1 = &input_tree;
// et T2 recoit l'arbre de reference
  T2 = &reference_tree;
  int r_size = T1->getDepth()*T1->getDegree();
  int c_size = T2->getNbVertex();
  int s_size = T1->getNbVertex();
  _treeConnectedTable.resize(1+r_size, 1+c_size, 1+s_size);
  _forestConnectedTable.resize(r_size, c_size, s_size);
  if (!T2->isNull()) 
  {
    _treeConnectedTable.openIntVector(T1->getNbVertex());
  }
}

// ------------------------------------------------
// Renvoie le nbre de cc dun arbre de racine
// input_vertex 
// ------------------------------------------------
int MatchingConnectedTable::getNBC(int input_vertex,int reference_vertex) const 
{
  int nbConnectedComponents;
// Si le noeud consideres est vide alors le nbre de cc est nul
  if ((input_vertex==EMPTY_NODE)||(reference_vertex==EMPTY_NODE)) {nbConnectedComponents=0; }
  else
    {
      nbConnectedComponents=_treeConnectedTable.getInt(input_vertex,reference_vertex);
    }
// on renvoie la  nbre  trouvee
  return(nbConnectedComponents);
}

int MatchingConnectedTable::getNBCF(int input_vertex,int reference_vertex) const 
{
  int nbConnectedComponents;
// Si le noeud consideres est vide alors le nbre de cc est nul
  if ((input_vertex==EMPTY_NODE)||(reference_vertex==EMPTY_NODE)) {nbConnectedComponents=0; }
  else
    {
      nbConnectedComponents=_forestConnectedTable.getInt(input_vertex,reference_vertex);
    }
// on renvoie la  nbre  trouvee
  return(nbConnectedComponents);
}


// ----------------------------------------
// Insere le nbre de cc entre deux arbres
// dans le tableau des cc entre arbres
// -----------------------------------------
void  MatchingConnectedTable::putNBC(int input_vertex,int reference_vertex ,int new_int)
{
  _treeConnectedTable.putInt(new_int,input_vertex,reference_vertex);
}

void  MatchingConnectedTable::putNBCF(int input_vertex,int reference_vertex ,int new_int)
{
  _forestConnectedTable.putInt(new_int,input_vertex,reference_vertex);
}

void MatchingConnectedTable::openIntVector(int input_vertex)
{
  _treeConnectedTable.openIntVector(input_vertex);
  _treeConnectedTable.putInt(0,input_vertex,T2->getNbVertex());
  _forestConnectedTable.openIntVector(input_vertex);
}


void MatchingConnectedTable::closeIntVector(int input_vertex)
{
  _treeConnectedTable.closeIntVector(input_vertex);
  _forestConnectedTable.closeIntVector(input_vertex);
}


