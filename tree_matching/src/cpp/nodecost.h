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


#ifndef SB_NODE_COST_HEADER
#define SB_NODE_COST_HEADER

#include <string>
using std::string;
#include <iostream>
#include <map>
using namespace std;

#include "definitions.h"
#include "treenode.h"

/** Type of  Node Cost */
enum NodeCostType { MATRIX, WEIGTH , TOPOLOGIC, LOCAL_TOPO, LOCAL_WEIGHT , SCORE};

/** Two type of Norms for computing the Node Cost */
enum Norm {L1, L2};


/**
 *\class NodeCost
 *\brief Definition of an elementary cost for the comparison between nodes
 *\author Pascal ferraro
 *\date 1999
 */

class NodeCost
{
  public :

    /** Default constructor. */
    NodeCost(){}

    /** Constructs a NodeCost with the type /e type and the default norm /e L1. */
    NodeCost(NodeCostType );

    NodeCost(char* );

    /** Constructs a NodeCost with the type /e type and the default norm /e norm . */
    NodeCost(NodeCostType,Norm );

    /** Destructor. */
    ~NodeCost(){}

    /** Returns the insertion cost of /e node */
    virtual DistanceType getInsertionCost(TreeNode* );

    int getCost(const char*, const char* );

    /** Returns the deletion cost of /e node */
    virtual DistanceType getDeletionCost(TreeNode* );

    /** Returns the changing cost between /e i_node and /e r_node*/
    virtual DistanceType getChangingCost(TreeNode* ,TreeNode* );

    /** Returns the type of the node cost*/
    NodeCostType type() const { return(_type); }

protected :
  
  NodeCostType _type;
  Norm         _norm;
  
  std::map<string, int> symbols;  
  int matrix [30][30];
  int matrix_size;
    
};

#endif

