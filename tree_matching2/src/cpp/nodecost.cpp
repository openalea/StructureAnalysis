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

/* ----------------------------------------------------------------------- */
NodeCost::~NodeCost(){ }

DistanceType NodeCost::getInsertionCost(const TreeNodePtr node) const
{ return 1; }

DistanceType NodeCost::getDeletionCost(const TreeNodePtr node) const
{ return 1; }

DistanceType NodeCost::getChangingCost(const TreeNodePtr i_node,const TreeNodePtr r_node) const
{ return 0; }

DistanceType NodeCost::getMergingCost(const vector<TreeNodePtr> i_node,const TreeNodePtr r_node) const
{ return 1; }

DistanceType NodeCost::getSplittingCost(const TreeNodePtr  i_node,const vector<TreeNodePtr> r_node) const
{ return 1; }
/* ----------------------------------------------------------------------- */


ScoreNodeCost::ScoreNodeCost( ):
NodeCost() { }

DistanceType ScoreNodeCost::getInsertionCost(const TreeNodePtr node) const
{ return -1; }

DistanceType ScoreNodeCost::getDeletionCost(const TreeNodePtr node) const
{ return -1; }

DistanceType ScoreNodeCost::getChangingCost(const TreeNodePtr i_node,const TreeNodePtr r_node) const
{ return 2; }

DistanceType ScoreNodeCost::getMergingCost(const vector<TreeNodePtr> i_node,const TreeNodePtr r_node) const
{ return -1; }

DistanceType ScoreNodeCost::getSplittingCost(const TreeNodePtr  i_node,const vector<TreeNodePtr> r_node) const
{ return -1; }

/* ----------------------------------------------------------------------- */
