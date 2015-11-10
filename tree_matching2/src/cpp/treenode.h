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
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>

#include "definitions.h"
using namespace std;


class TreeNode;
// A smart pointer on a TreeNode
typedef boost::shared_ptr<TreeNode> TreeNodePtr;
// A weak reference to avoid loop.
typedef boost::weak_ptr<TreeNode> TreeNodeWeakPtr;

#define NOID -1

/**
 *\class TreeNode
 *\brief Definition of a vertex in a tree graph
 *\author Pascal ferraro
 *\date 2009
 */


class TREEMATCH_API TreeNode
{

  public :
	/// A type for attributes of TreeNodes.
	typedef boost::any ValueType;
	typedef std::vector<ValueType> ValueVector;	

	typedef int IdType;
	typedef std::vector<IdType> ChildVector;

	// static const int NOID;

    /** Default constructor
     *\par Remarks
     * Never used
     */
    TreeNode() : _id(0), _father(NOID), _values() {};



  /** constructs a TreeNode with  /e id, /e father */
  TreeNode(IdType ,IdType father = NOID);

  /** Destructor. */
  ~TreeNode();

  /** \par Reading Functions. */

  /** Returns the father of TreeNode. */
  inline attribute_deprecated IdType father() const { return(_father);}
  inline IdType getFather() const { return _father;}

  /** Returns the id of the TreeNode in MTG.
   *  id is the reference of TreeNode in a TreeGraph */
  inline IdType getId() const { return _id;}

    /** \par Writing Functions. */

  /** Set the id of the TreeNode in a TreeGraph.
   *  id is the reference of TreeNode in a TreeGraph */
  inline attribute_deprecated void putId(IdType id) { _id=id;}
  inline void setId(IdType id) { _id=id;}
   
  /** Set the father of the TreeNode in MTG. */
  inline attribute_deprecated void putFather(IdType father ){ _father=father;}
  inline void setFather(IdType father ){ _father=father;}
  
  /** Print a TreeNode*/
  void print() const;
  
  /** Returns the number of values into the TreeNode.*/
  inline size_t getValueSize() const {  return _values.size(); }

  /** Returns the values into the TreeNode.*/
  inline const ValueVector& getValueList() const {  return _values; }
  inline ValueVector& getValueList() {  return _values; }

  /// Return untyped value 
  inline boost::any getAnyValue(size_t index) const 
  { assert(index<_values.size()); return _values[index]; }

  /// check if value at index is of type T
  template <class T>
  inline bool isValueOfType(size_t index) const 
  { assert(index<_values.size()); return _values[index].type() == typeid(T); }

  /// Return value of a given type
  template <class T>
  T getTypedValue(size_t index) const 
  { assert(index<_values.size()); return boost::any_cast<T>(_values[index]); }

  /// Set a value of any type at index
  template <class T>
  void setTypedValue(size_t index, const T& value) 
  { assert(index<_values.size()); _values[index]=boost::any(value); }

  /// Append a value of any type
  template <class T>
  void appendTypedValue(const T& value) 
  { _values.push_back(boost::any(value)); }

  DistanceType getValue(int index) const;

  /** Puts a value /e new_value at the index /e index to a TreeNode.*/
  inline attribute_deprecated void putValue(size_t index, DistanceType new_value = 0.) { setTypedValue(index,new_value); }

 /** Puts a value /e new_value at the index /e index to a TreeNode.*/
  inline attribute_deprecated void addValue( DistanceType new_value = 0.) { appendTypedValue(new_value); }

  /** Add a child in child list, a child is referenced by is id only */
  inline void addChild(IdType child_id){ _childList.push_back(child_id); }

  /** Set Child List */
  inline void setChildList(const ChildVector& child_list) { _childList = child_list; }

  /** Get Child List */
  inline const ChildVector& getChildList() const { return _childList; }

  /** Get i-th Child  */
  IdType getChild(int i) const;

  /** Get Number of Child */
  inline size_t  getChildNumber() const { return _childList.size(); }
  inline bool  hasChild() const { return !_childList.empty(); }

  /** depth accessors */
  inline attribute_deprecated int depth() const{ return _depth;  }  
  inline attribute_deprecated void putDepth(int depth){ _depth = depth;  }

  inline int getDepth() const{ return _depth;  }  
  inline void setDepth(int depth){ _depth = depth;  }
  
  public:
  /** 
    This class make it possible to construct TreeNode.
	The interesting point is that it allow to redefine its behaviour
	by passing a custom builder function. By default, it uses the TreeNode constructor.
	Note that it is singleton.
    */
  class TREEMATCH_API Factory {
  public:
	  typedef TreeNodePtr (*TNBuilder)(IdType, IdType);
	  friend class TreeNode;

	  // set of a custion treenode builder function 
	  inline void setBuilder(TNBuilder b) { __builder = b; }
	  // reset builder function
	  inline void setBuilderToDefault() { __builder = NULL; }

	  TreeNodePtr build(IdType, IdType);
  protected:
	   Factory();
	   TNBuilder __builder;
  };

  // singleton access
  static Factory& factory();

  private :
  IdType _id;
  IdType _father;
  int _depth;
  ChildVector _childList;
  ValueVector _values;
};


#endif

