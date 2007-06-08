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


#include "emdtable.h"



ExtendedMatchingDistanceTable::ExtendedMatchingDistanceTable(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _inputTreeToEmpty.resize(input.getNbVertex());
  _referenceTreeFromEmpty.resize(reference.getNbVertex());
  int data_size=input.getDepth()*input.getDegree();
  _treeDistTable.resize(1+data_size,1+reference.getNbVertex(),1+input.getNbVertex());
  _forestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _minusForestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _plusForestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _imfrfDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _ifrmfDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  if (!T1->isNull()) inputTreeToEmpty(0);
  if (!T2->isNull()) 
  {
    _treeDistTable.openDistanceVector(T1->getNbVertex());
    referenceTreeFromEmpty(0);
  }
}


void ExtendedMatchingDistanceTable::make(TreeGraph& input,TreeGraph& reference,NodeCost& nodeDistance)
{
  T1=&input;
  T2=&reference;
  ND=&nodeDistance;
  _inputTreeToEmpty.resize(input.getNbVertex());
  _referenceTreeFromEmpty.resize(reference.getNbVertex());
  int data_size=input.getDepth()*input.getDegree();
  _treeDistTable.resize(1+data_size,1+reference.getNbVertex(),1+input.getNbVertex());
  _forestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _minusForestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _plusForestDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _imfrfDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  _ifrmfDistTable.resize(data_size,reference.getNbVertex(),input.getNbVertex());
  if (!T1->isNull()) inputTreeToEmpty(0);
  if (!T2->isNull()) 
  {
    _treeDistTable.openDistanceVector(T1->getNbVertex());
    referenceTreeFromEmpty(0);
  }
}

DistanceType ExtendedMatchingDistanceTable::inputTreeToEmpty(int vertex)
{
  DistanceType cost=getDCost(vertex);
  if (!T1->isLeaf(vertex))
  {
    cost=cost+inputForestToEmpty(vertex);
  }
  assert(cost>=0);
  _inputTreeToEmpty[vertex]=cost;
  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::inputForestToEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=1;i<=T1->getNbChild(vertex);i++)
  {
    cost=cost+inputTreeToEmpty(T1->child(vertex,i));
  }
  assert(cost>=0);
  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::referenceTreeFromEmpty(int vertex)
{
  DistanceType cost=getICost(vertex);
  if (!T2->isLeaf(vertex))
  {
    cost=cost+referenceForestFromEmpty(vertex);
  }
  assert(cost>=0);
  _referenceTreeFromEmpty[vertex]=cost;
  _treeDistTable.putDistance(cost,T1->getNbVertex(),vertex);
  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::referenceForestFromEmpty(int vertex)
{
  DistanceType cost=MINDIST;
  for (int i=1;i<=T2->getNbChild(vertex);i++)
  {
    cost=cost+referenceTreeFromEmpty(T2->child(vertex,i));
  }
  assert(cost>=0);
  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::getDBT(int input_vertex,int reference_vertex) const 
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
    distance=_referenceTreeFromEmpty[reference_vertex];
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      distance=_inputTreeToEmpty[input_vertex];
    }
    else
    {
      distance=_treeDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

DistanceType ExtendedMatchingDistanceTable::getDBF(int input_vertex,int reference_vertex) const 
{
  DistanceType distance;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
    distance=_referenceTreeFromEmpty[reference_vertex]-getICost(reference_vertex);
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      distance=_inputTreeToEmpty[input_vertex]-getDCost(input_vertex);
    }
    else
    {
      distance=_forestDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

DistanceType ExtendedMatchingDistanceTable::getDBPF(int input_vertex,int reference_vertex) const 
{
  DistanceType distance=MINDIST;
  int son;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
    for (son=1;son<=T2->getNbChild(reference_vertex);son++) 
    { 
      int reference_son=T2->child(reference_vertex,son);
      if (sameComplex(T2->getNode(reference_son),T2->getNode(reference_vertex)))
      {
	distance=distance+_referenceTreeFromEmpty[reference_son];
      }
    }
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      for (son=1;son<=T1->getNbChild(input_vertex);son++) 
      { 
	int input_son=T1->child(input_vertex,son);
	if (sameComplex(T1->getNode(input_son),T1->getNode(input_vertex)))
	{
	  distance=distance+_inputTreeToEmpty[input_son];
	}
      }
    }
    else
    {
      distance=_plusForestDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

DistanceType ExtendedMatchingDistanceTable::getDBMF(int input_vertex,int reference_vertex) const 
{
  DistanceType distance=MINDIST;
  int son;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
    for (son=1;son<=T2->getNbChild(reference_vertex);son++) 
    { 
      int reference_son=T2->child(reference_vertex,son);
      if (!sameComplex(T2->getNode(reference_son),T2->getNode(reference_vertex)))
      {
	distance=distance+_referenceTreeFromEmpty[reference_son];
      }
    }
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      for (son=1;son<=T1->getNbChild(input_vertex);son++) 
      { 
	int input_son=T1->child(input_vertex,son);
	if (!sameComplex(T1->getNode(input_son),T1->getNode(input_vertex)))
	{
	  distance=distance+_inputTreeToEmpty[input_son];
	}
      }
    }
    else
    {
      distance=_minusForestDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

DistanceType ExtendedMatchingDistanceTable::getDBIMFRF(int input_vertex,int reference_vertex) const 
{
  DistanceType distance=MINDIST;
  int son;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
     distance=_referenceTreeFromEmpty[reference_vertex]-getICost(reference_vertex);
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      for (son=1;son<=T1->getNbChild(input_vertex);son++) 
      { 
	int input_son=T1->child(input_vertex,son);
	if (!sameComplex(T1->getNode(input_son),T1->getNode(input_vertex)))
	{
	  distance=distance+_inputTreeToEmpty[input_son];
	}
      }
    }
    else
    {
      distance=_imfrfDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

DistanceType ExtendedMatchingDistanceTable::getDBIFRMF(int input_vertex,int reference_vertex) const 
{
  DistanceType distance=MINDIST;
  int son;
  if ((input_vertex==EMPTY_NODE)&&(reference_vertex==EMPTY_NODE)) {distance=MINDIST;};
  if (input_vertex==EMPTY_NODE)
  {
    for (son=1;son<=T2->getNbChild(reference_vertex);son++) 
    { 
      int reference_son=T2->child(reference_vertex,son);
      if (!sameComplex(T2->getNode(reference_son),T2->getNode(reference_vertex)))
      {
	distance=distance+_referenceTreeFromEmpty[reference_son];
      }
    }
  }
  else
  {
    if (reference_vertex==EMPTY_NODE)
    {
      distance=_referenceTreeFromEmpty[input_vertex]-getICost(input_vertex);
    }
    else
    {
      distance=_ifrmfDistTable.getDistance(input_vertex,reference_vertex);
    }
  }
  return(distance);
}

void  ExtendedMatchingDistanceTable::putDBF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_forestDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void  ExtendedMatchingDistanceTable::putDBT(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_treeDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void  ExtendedMatchingDistanceTable::putDBMF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_minusForestDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void  ExtendedMatchingDistanceTable::putDBPF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_plusForestDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void  ExtendedMatchingDistanceTable::putDBIMFRF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_imfrfDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void  ExtendedMatchingDistanceTable::putDBIFRMF(int input_vertex,int reference_vertex ,DistanceType new_distance)
{
	_ifrmfDistTable.putDistance(new_distance,input_vertex,reference_vertex);
}

void ExtendedMatchingDistanceTable::openDistancesVector(int input_vertex)
{
  _treeDistTable.openDistanceVector(input_vertex);
  _treeDistTable.putDistance(_inputTreeToEmpty[input_vertex],input_vertex,T2->getNbVertex());
  _forestDistTable.openDistanceVector(input_vertex);
  _plusForestDistTable.openDistanceVector(input_vertex);
  _minusForestDistTable.openDistanceVector(input_vertex);
  _ifrmfDistTable.openDistanceVector(input_vertex);
  _imfrfDistTable.openDistanceVector(input_vertex);
}


void ExtendedMatchingDistanceTable::closeDistancesVector(int input_vertex)
{
  _treeDistTable.closeDistanceVector(input_vertex);
  _forestDistTable.closeDistanceVector(input_vertex);
  _plusForestDistTable.closeDistanceVector(input_vertex);
  _minusForestDistTable.closeDistanceVector(input_vertex);
  _ifrmfDistTable.closeDistanceVector(input_vertex);
  _imfrfDistTable.closeDistanceVector(input_vertex);
}

DistanceType ExtendedMatchingDistanceTable::getICost(int vertex) const
{
  TreeNode* node=T2->getNode(vertex);
  DistanceType cost;

  switch(ND->type())
  { 
      case TOPOLOGIC  : cost=ND->getInsertionCost(node); break;
      case WEIGTH     : cost=((WeightedNodeCost*) ND)->getInsertionCost(node);break;
      case MATRIX     : cost=((MatrixNodeCost*)   ND)->getInsertionCost(node);break;
      default         : assert(0);break;
  }

  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::getDCost(int vertex ) const
{
  TreeNode* node=T1->getNode(vertex);
  DistanceType cost;
  
  switch(ND->type())
  { 
     case TOPOLOGIC : { cost = ND->getDeletionCost(node); };break;
     case WEIGTH    : { cost = ((WeightedNodeCost*) ND)->getDeletionCost(node); };break;
     case MATRIX    : { cost = ((MatrixNodeCost*) ND)->getDeletionCost(node); };break;
     default        : { assert(0); };break;
  }
  
  return(cost);
}

DistanceType ExtendedMatchingDistanceTable::getCCost(int vertex1 ,int vertex2) const
{
	TreeNode* node1=T1->getNode(vertex1);
	TreeNode* node2=T2->getNode(vertex2);
	DistanceType cost;

	switch(ND->type())
       	{ 
        	case TOPOLOGIC  : {cost=ND->getChangingCost(node1,node2);};break;
        	case WEIGTH     : {cost=((WeightedNodeCost*) ND)->getChangingCost(node1,node2);};break;
        	case MATRIX     : {cost=((MatrixNodeCost*) ND)->getChangingCost(node1,node2);};break;
        	default         : assert(0);break;
        }
	
	return(cost);
}

Boolean ExtendedMatchingDistanceTable::sameComplex(TreeNode* tree_node1,TreeNode* tree_node2) const
{
  if ( tree_node1->getComplex() == tree_node2->getComplex())
  {
    return(TRUE);
  }
  else
  {
    return(FALSE);
  };
}
















