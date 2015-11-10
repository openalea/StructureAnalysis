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


#include"wnodecostnonconst.h"


WeightedNodeCostNonConst::WeightedNodeCostNonConst(NodeCostType type ,WeightList weights,DistanceType ins_del_cost,Norm norm)
{
  _type=type;
  _weights=weights;
  _InsDelCost=ins_del_cost;
  _norm=norm;
}

WeightedNodeCostNonConst::WeightedNodeCostNonConst(NodeCostType type ,WeightList weights,DistanceType ins_del_cost)
{
  _type=type;
  _weights=weights;
  _InsDelCost=ins_del_cost;
  _norm=L1;
}

DistanceType WeightedNodeCostNonConst::getInsertionCost(TreeNode* node)
{
  DistanceType cost=0;
  VId vertex1=node->getVertex();
  if (_norm==L1)
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if (f_value1==DIST_UNDEF)f_value1=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*ABS(f_value1);
	}
    }
  else
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*f_value1*f_value1;
	}
    }
  assert(cost>=0);
  cost=sqrt(cost);
  return(cost);
}

DistanceType WeightedNodeCostNonConst::getDeletionCost(TreeNode* node)
{
  DistanceType cost=0;
  VId vertex1=node->getVertex();
  if (_norm==L1)
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*ABS(f_value1);
	}
    }
  else
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*f_value1*f_value1;
	}
    }
  
  assert(cost>=0);
  cost=sqrt(cost);
  return(cost);
}

DistanceType WeightedNodeCostNonConst::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost=0;
  VId vertex1=i_node->getVertex();
  VId vertex2=r_node->getVertex();
  if (_norm==L1)
    {
      for (int i=0;i<i_node->getValueSize();i++) 
	{
	  DistanceType f_value1=i_node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  DistanceType f_value2=r_node->getValue(i);
	  if (f_value2==DIST_UNDEF) f_value2=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*ABS(f_value1-f_value2);
	}
    }
  else
    {
      for (int i=0;i<i_node->getValueSize();i++) 
	{
	  DistanceType f_value1=i_node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  DistanceType f_value2=r_node->getValue(i);
	  if (f_value2==DIST_UNDEF) f_value2=0;
	  list<DistanceType>::iterator begin;
	  begin = _weights.begin();
	  for (int j=0;j<i;j++)
	    begin++;
	  cost=cost+(*begin)*(f_value1-f_value2)*(f_value1-f_value2);
	}
    }
  assert(cost>=0);
  cost=sqrt(cost);
  return(cost);
}
 












