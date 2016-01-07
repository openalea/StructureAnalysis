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
 *       $Id: nodecost.h 3258 2007-06-06 13:18:26Z dufourko $
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


#include "definitions.h"
#include "treenode.h"
#include <boost/shared_ptr.hpp>

/* ----------------------------------------------------------------------- */


/**
 *\class NodeCost
 *\brief Definition of an elementary cost for the comparison between nodes
 *\author Pascal ferraro
 *\date 2009
 */

class TREEMATCH_API NodeCost
{
public :

  /** Default constructor. */
  NodeCost(){}

  /** Destructor. */
  virtual ~NodeCost();
  
  /** Returns the insertion cost of /e node */
  virtual DistanceType getInsertionCost(const TreeNodePtr ) const ;
  
  /** Returns the deletion cost of /e node */
  virtual DistanceType getDeletionCost( const TreeNodePtr ) const ;
  
  /** Returns the changing cost between /e i_node and /e r_node*/
  virtual DistanceType getChangingCost( const TreeNodePtr , const TreeNodePtr ) const ;

  /** Returns the merging cost between /e a vector of node and /e r_node*/
  virtual DistanceType getMergingCost( const vector<TreeNodePtr> , const TreeNodePtr ) const ;
       
  /** Returns the splitting cost between /e a vector of node and /e r_node*/
  virtual DistanceType getSplittingCost( const TreeNodePtr, const vector<TreeNodePtr> ) const ;
       
};


// A smart pointer on a NodeCost
typedef boost::shared_ptr<NodeCost> NodeCostPtr;
// A weak reference to avoid loop.
typedef boost::weak_ptr<NodeCost> NodeCostWeakPtr;


/* ----------------------------------------------------------------------- */


class TREEMATCH_API ScoreNodeCost : public NodeCost
{
public :

    /** Constructs a NodeCost with the type /e type and the default norm /e norm . */
  ScoreNodeCost();

  /** Destructor. */
  virtual ~ScoreNodeCost(){}
  
  /** Returns the insertion cost of /e node */
  virtual DistanceType getInsertionCost(const TreeNodePtr ) const ;
  
  /** Returns the deletion cost of /e node */
  virtual DistanceType getDeletionCost(const TreeNodePtr ) const ;
  
  /** Returns the changing cost between /e i_node and /e r_node*/
  virtual DistanceType getChangingCost(const TreeNodePtr ,const TreeNodePtr ) const ;

  /** Returns the merging cost between /e a vector of node and /e r_node*/
  virtual DistanceType getMergingCost( const vector<TreeNodePtr> , const TreeNodePtr ) const ;
         
  /** Returns the splitting cost between /e a vector of node and /e r_node*/
  virtual DistanceType getSplittingCost( const TreeNodePtr, const vector<TreeNodePtr> ) const ;
};

/* ----------------------------------------------------------------------- */


#endif

