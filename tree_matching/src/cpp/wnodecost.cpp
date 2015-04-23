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


#include"wnodecost.h"

using namespace stat_tool;


WeightedNodeCost::WeightedNodeCost(NodeCostType type, const VectorDistance& ivect, const ValueVector& dispersion,
                                   const ValueVector& max, const ValueVector& min, DistanceType indelcoeff)

{
  _type=type;
  _vectDist=*(new VectorDistance(ivect));
  _dispersion=dispersion;
  _maxValue=max;
  _minValue=min;
  _InsDelCostCoeff=indelcoeff;
}


DistanceType WeightedNodeCost::getInsertionCost(TreeNode* node)
{
DistanceType res =0;
  DistanceType cost=0;
  if (_vectDist.get_distance_type()==ABSOLUTE_VALUE)
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if  (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])/_dispersion[i];
	    }
	  else
	    cost=cost+_vectDist.get_weight(i)*_vectDist.max_symbol_distance_computation(i)[int(f_value1)]/_dispersion[i];
	}
      assert(cost>=0);
    }
  else
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if  (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])
		*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])/_dispersion[i];
	    }
	  else
	    cost=cost+_vectDist.get_weight(i)*_vectDist.max_symbol_distance_computation(i)[int(f_value1)]
	      *_vectDist.max_symbol_distance_computation(i)[int(f_value1)]/_dispersion[i];
	}
      assert(cost>=0);
      cost=sqrt(cost);
    }

  if(_type==SCORE)
  {res =-(_InsDelCostCoeff*cost);
  }
else
  res=_InsDelCostCoeff*cost;
return res;
}

DistanceType WeightedNodeCost::getDeletionCost(TreeNode* node)
{
DistanceType res = 0;
  DistanceType cost=0;
  if (_vectDist.get_distance_type()==ABSOLUTE_VALUE)
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if  (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])/_dispersion[i];
	    }
	  else
	    cost=cost+_vectDist.get_weight(i)*_vectDist.max_symbol_distance_computation(i)[int(f_value1)]/_dispersion[i];
	}
      assert(cost>=0);
    }
  else
    {
      for (int i=0;i<node->getValueSize();i++) 
	{
	  DistanceType f_value1=node->getValue(i);
	  if  (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])
		*MAX(_maxValue[i]-f_value1,f_value1-_minValue[i])/_dispersion[i];
	    }
	  else
	    cost=cost+_vectDist.get_weight(i)*_vectDist.max_symbol_distance_computation(i)[int(f_value1)]
	      *_vectDist.max_symbol_distance_computation(i)[int(f_value1)]/_dispersion[i];
	}
      assert(cost>=0);
      cost=sqrt(cost);
    }
if(_type==SCORE)
{res =-(_InsDelCostCoeff*cost);
}
else
  res=_InsDelCostCoeff*cost;
return res;

}

DistanceType WeightedNodeCost::getChangingCost(TreeNode* i_node,TreeNode* r_node)
{
DistanceType res = 0;
  DistanceType cost=0;
  if (_vectDist.get_distance_type()==ABSOLUTE_VALUE)
    {
      for (int i=0;i<i_node->getValueSize();i++) 
	{
	  DistanceType f_value1=i_node->getValue(i);
	  if (f_value1==DIST_UNDEF) f_value1=0;
	  DistanceType f_value2=r_node->getValue(i);
	  if (f_value2==DIST_UNDEF) f_value2=0;

	  if (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*ABS(f_value1-f_value2)/_dispersion[i];
	      //cout<<ABS(f_value1-f_value2)/_dispersion[i]<<" - "<<f_value1<<" - "<<f_value2<<" - "<<_dispersion[i]<<endl;
	    }
	  else
	    {
	      cost=cost+_vectDist.get_weight(i)
		*_vectDist.get_symbol_distance(i,int(f_value1),int(f_value2))/_dispersion[i];
	    }
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

	  if (_vectDist.get_variable_type(i)==NUMERIC)
	    {
	      cost=cost+_vectDist.get_weight(i)*(f_value1-f_value2)*(f_value1-f_value2)/_dispersion[i];
	    }
	  else
 	      cost=cost+_vectDist.get_weight(i)
		*_vectDist.get_symbol_distance(i,int(f_value1),int(f_value2))
		*_vectDist.get_symbol_distance(i,int(f_value1),int(f_value2))
		/_dispersion[i];		       
	}
      cost=sqrt(cost);
    }
  assert(cost>=0);
if(_type==SCORE)
{res = 2-cost;
}
else{
  res=cost;
}
return res;
}
 












