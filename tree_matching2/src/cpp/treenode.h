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
 *       $Id: treenode.h 3264 2007-06-06 14:22:22Z dufourko $
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


#ifndef SB_TREE_NODE_HEADER
#define SB_TREE_NODE_HEADER

#include <list>
#include <iostream>
#include "definitions.h"
using namespace std;
typedef DistanceType ValueType;
/**
 * A type for attributes of TreeNodes.
 */
typedef std::vector<DistanceType> ValueVector;


/**
 *\class TreeNode
 *\brief Definition of a vertex in a tree graph
 *\author Pascal ferraro
 *\date 2009
 */


class TreeNode
{

  public :

    /** Default constructor
     *\par Remarks
     * Never used
     */
    TreeNode() : _id(0),
		 _father(0),
		 _values(NULL) {};

//     /** Copy constructor. */
//   TreeNode(const TreeNode& inode):  
//     _id(inode.getId()),
//     _father(inode.father()),
//     _values(NULL){
//     //    cout <<"_nbvalue =" <<_nbvalue<<endl;
//     for(int i=0;i<_values.size();i++)
//       putValue(i,inode.getValue(i));
//   }


  /** constructs a TreeNode with  /e id, /e father */
  TreeNode(int ,int);

  /** Destructor. */
  ~TreeNode();

  /** \par Reading Functions. */

  /** Returns the father of TreeNode. */
  int father() const { return(_father);}

  /** Returns the id of the TreeNode in MTG.
   *  id is the reference of TreeNode in a TreeGraph */
  int getId() const { return(_id);}

    /** \par Writing Functions. */

  /** Puts the id of the TreeNode in a TreeGraph.
   *  id is the reference of TreeNode in a TreeGraph */
  void putId(int id) { _id=id;}
   
  /** Puts the father of the TreeNode in MTG. */
  void putFather(int father ){ _father=father;}
  
  /** Print a TreeNode*/
  void print() const;
  
  /** Returns the id of value affected to a TreeNode.*/
  int getValueSize() const;
  DistanceType getValue(int index) const;

  /** Puts a value /e new_value at teh index /e index to a TreeNode.*/
  void putValue(int index, DistanceType new_value = 0.);
 /** Puts a value /e new_value at teh index /e index to a TreeNode.*/
  void addValue( DistanceType new_value = 0.);

  /** Add a child in child list, a child is referenced by is id only */
  void addChild(int);

  /** Set Child List */
  void setChildList(vector<int>);

  /** Get Child List */
  vector<int> getChildList() const;

  /** Get Child  */
  int getChild(int) const;

  /** Get Child Number */
  int getChildNumber() const;
  int depth() const{
    return _depth;
  }
  
  void putDepth(int depth){
    _depth = depth;
  }
  
  private :
  int _id;
  int _father;
  int _depth;
  vector<int> _childList;
  ValueVector _values;
};

#endif

