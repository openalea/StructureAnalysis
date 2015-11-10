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


#include"melodyNodeCost.h"


MelodyNodeCost::MelodyNodeCost(NodeCostType type,ValueVector weights)

{
  _type=type;
  _InsDelCost = 0.19;
  _cost_pitch.resize(weights.size());
  for (int i =0;i<8;i++){
  _cost_pitch[i] = weights[i];
  }	
  /*_cost_pitch[0] = 0:
  _cost_pitch[1] = 1.9*4;
  _cost_pitch[2] = 1.65*4;
  _cost_pitch[3] = 0.55;
  _cost_pitch[4] = 0.55;
  _cost_pitch[5] = 0.2;
  _cost_pitch[6] = 1.2*4;
  _cost_pitch[7] = 0.1;*/
}

int modulo12(int n)
{
  int a = n%12;

  if (a < -6)
    a+=12;
  if (a >= 6)
    a-=12;

  return a;
}

DistanceType MelodyNodeCost::getInsertionCost(TreeNode* node)
{  
	DistanceType cost=0.;
  //prise en compte du  time
  DistanceType f_value=node->getValue(0);
  cost=cost+ABS(f_value)* _InsDelCost;
//  cerr<<node->getVertex()<<" - "<<cost<<endl;
  if (f_value=node->getValue(1)!=8)
      cost += _cost_pitch[7];

  assert(cost>=0);
  return(cost);

}

DistanceType MelodyNodeCost::getDeletionCost(TreeNode* node)
{
  return(getInsertionCost(node));
}

DistanceType MelodyNodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
  DistanceType cost=0;
  //prise en compte du  time
  DistanceType f_value1=i_node->getValue(0);
  DistanceType f_value2=r_node->getValue(0);
  cost=cost+ABS(f_value1-f_value2)* _InsDelCost;

  //prise en compte du pitch
  
  f_value1=i_node->getValue(1);
  f_value2=r_node->getValue(1);

  int pitch_variation=abs(modulo12(int(f_value1) - int(f_value2)));

  if (pitch_variation<=6)
    cost += _cost_pitch[pitch_variation];
  else 
    cost += _cost_pitch[7];
  assert(cost>=0);
  return(cost);
}
 











