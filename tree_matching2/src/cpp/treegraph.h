/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 *
 *       $Source: /usr/cvsmaster/AMAPmod/src/TreeMatching/treegraph.h,v $
 *       $Id: treegraph.h 3264 2007-06-06 14:22:22Z dufourko $
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


#ifndef SB_TREE_GRAPH_HEADER
#define SB_TREE_GRAPH_HEADER

#include "definitions.h"
#include "treenode.h"
#include <list>
#include <vector>
#include <iterator>


/**
 * List of TreeNode
 */
typedef std::vector<TreeNode> TreeNodeList;

/**
 * List of Children
 */

typedef std::vector<int> NodeList;


using namespace std;


/**
 *\class TreeGraph
 *\brief Definition of a tree graph
 *\author Pascal ferraro
 *\date 2009
 */

class TREEMATCH_API TreeGraph
{

  public :
  typedef TreeNodeList::const_iterator const_iterator;
    /** Constructor */
    TreeGraph();
    /** Destructor */
    ~TreeGraph();

  void addNode(int, int);

  void addNode(TreeNode);

  /** return TreeNode corresponding to vertex */
  TreeNode* getTreeNode(int vertex) ;

   /** Return the father of a node */
  int father(int ) const;

  /** Give a child of the ChildrenList */
  int child(const int ,const int) const ;

  void addValue(int, DistanceType);
  DistanceType getValue(int, int) const;


  /** Give the SonsList */
  const NodeList childList(const int ) const ;

  /** Return the number of child of a node  */
  int getNbChild(int ) const;
  int getNbDesc(int ) const;

  /** Return the number of node in the tree */
  int getNbVertex() const;

  /** Return the root id */
  int getRoot() const{
    return _rootId;
  }

  /** Return the root id */
  int getDegree() const{
    return _degree;
  }
  
  /** Return whether the node is a leaf or not          */
  int isLeaf(int ) const;

  /** Return a node */
  TreeNode& getNode(int );

  int getNumber(int vertex) const;


  /** Return wheter the tree is null or not  */
  int isNull();

  /** Method for printing */
  void print() const;
  ostream& mtg_write(ostream &os) const;
  bool mtg_write(  char *path ) ;   



  private :
  int _nbVertex;
  TreeNodeList _treenodes;
  int _depth;
  int _degree;
  int _rootId;
};


#endif


