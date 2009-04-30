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
 *       $Id: nodecost.cpp 3258 2007-06-06 13:18:26Z dufourko $
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


#include"nodecost.h"

NodeCost::NodeCost(string type){
  _type = TOPOLOGIC;
  _norm = L1;
}

NodeCost::NodeCost( NodeCostType type){
  _type=type;
  _norm = L1;
}

NodeCost::NodeCost( NodeCostType type,Norm norm){
  _type=type;
  _norm = norm;
}

DistanceType NodeCost::getInsertionCost(TreeNode* node){
  DistanceType cost;

  if (_type != SCORE)
    cost = 1;
  else
    cost= -1;
  return(cost);
}

DistanceType NodeCost::getDeletionCost(TreeNode* node)
{
  DistanceType cost;
   if (_type != SCORE)
     cost = 1;
   else
     cost= -1;
    return(cost);
}

DistanceType NodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost;
   if (_type != SCORE)
     cost = 0;
   else
     cost = 2;
  return(cost);
}


