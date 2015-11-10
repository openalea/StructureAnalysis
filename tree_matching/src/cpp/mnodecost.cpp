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


#include "mnodecost.h"

MatrixNodeCost::MatrixNodeCost(DistanceVectorTable& dtable,NodeCostType type)
{
  _dtable=&dtable;
  _type=type;
}

DistanceType MatrixNodeCost::getInsertionCost(TreeNode* node)
{
  DistanceType cost=getDist(-1,node->getNumber());
  return(cost);
}

DistanceType MatrixNodeCost::getDeletionCost(TreeNode* node)
{
  DistanceType cost=getDist(node->getNumber(),-1);
  return(cost);
}

DistanceType MatrixNodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost=getDist(i_node->getNumber(),r_node->getNumber());
  return(cost);
}

DistanceType MatrixNodeCost::getDist(int i_node,int r_node)
{
  if (i_node==-1) 
    {
      return((*_dtable)[(*_dtable).size()-1][r_node]);
    }
  else
    {
      if (r_node==-1) 
	{
	  return((*_dtable)[i_node][(*_dtable)[i_node].size()-1]);
	}
      else
	{
	  return((*_dtable)[i_node][r_node]);
	}
    }
}

