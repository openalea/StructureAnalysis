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


#ifndef SB_WEIGTHED_NODE_COST_HEADER
#define SB_WEIGTHED_NODE_COST_HEADER

#include "definitions.h"
#include "treenode.h"
#include "nodecost.h"
//#include "mdtable.h"
#include "stat_tool/stat_tools.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"

const DistanceType NOTDEF=-1.0;

typedef list<DistanceType> WeightList;

/**
 *\class WeightedNodeCost
 *\brief Definition of a weighted cost for the comparison between nodes
 * inherit class of NodeCost
 *\author Pascal ferraro
 *\date 1999
 */
class WeightedNodeCost : public NodeCost
{
  public :

    /** Default constructor. */
    WeightedNodeCost(){}

    /** Destructor. */
    ~WeightedNodeCost(){}

    /** Constructs a NodeCost with the type /e type and with a VectorDistance, a dispersion, a maximum and minimum values.and a product coefficient for indel cost. */
    WeightedNodeCost(NodeCostType,
		     const VectorDistance&,
		     const ValueVector&,
		     const ValueVector&,
		     const ValueVector&,DistanceType);

    /** Returns the insertion cost of /e node */    
    DistanceType getInsertionCost(TreeNode* );

    /** Returns the deletion cost of /e node */
    DistanceType getDeletionCost(TreeNode* );

    /** Returns the changing cost between /e i_node and /e r_node*/
    DistanceType getChangingCost(TreeNode* ,TreeNode* );

  private :
    WeightList       _weights;
    DistanceType     _InsDelCost;
    VectorDistance   _vectDist;
    ValueVector      _dispersion;
    ValueVector      _maxValue;
    ValueVector      _minValue;
    DistanceType     _InsDelCostCoeff;
};

#endif

