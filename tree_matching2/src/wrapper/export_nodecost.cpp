/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       TreeMatching : Comparison of Tree Structures
 *
 *       Copyright 1995-2009 UMR LaBRI
 *
 *       File author(s): P.ferraro (pascal.ferraro@labri.fr)
 * 
 # ---------------------------------------------------------------------------
 #
 #                      GNU General Public Licence
 #
 #       This program is free software; you can redistribute it and/or
 #       modify it under the terms of the GNU General Public License as
 #       published by the Free Software Foundation; either version 2 of
 #       the License, or (at your option) any later version.
 #
 #       This program is distributed in the hope that it will be useful,
 #       but WITHOUT ANY WARRANTY; without even the implied warranty of
 #       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 #       GNU General Public License for more details.
 #
 #       You should have received a copy of the GNU General Public
 #       License along with this program; see the file COPYING. If not,
 #       write to the Free Software Foundation, Inc., 59
 #       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 #
 # ---------------------------------------------------------------------------
 */


//#include "extract_list.h"
#include "export_list.h"
#include <boost/python.hpp>
#include "nodecost.h"

using namespace boost::python;
#define bp boost::python

#define treenode2python(node) bp::object(node)
inline object treenodelist2python( const vector<TreeNodePtr>& nodes )
{
  return make_list(nodes)();
 }

class PyNodeCost : public NodeCost, public bp::wrapper<NodeCost>
{
public:

	//Constructor
    PyNodeCost() :  NodeCost(), bp::wrapper<NodeCost>() {   }

   virtual ~PyNodeCost() {}

   /** Returns the insertion cost of /e node */
   virtual DistanceType getInsertionCost(const TreeNodePtr node) const
    { 
	  if(bp::override f = this->get_override("getInsertionCost"))
		return bp::call<DistanceType>(f.ptr(), treenode2python(node));
	  else return default_getInsertionCost(node);
	}

    /** Default implementation of insertion cost method. Required for boost python*/
    DistanceType default_getInsertionCost(const TreeNodePtr node) const
    {  return NodeCost::getInsertionCost(node); }

    /** Returns the deletion cost of /e node */
    virtual DistanceType getDeletionCost(const TreeNodePtr node) const
	{ 
	  if(bp::override f = this->get_override("getDeletionCost"))
		return bp::call<DistanceType>(f.ptr(), treenode2python(node));
	  else return default_getDeletionCost(node);
	}

    /** Default implementation of deletion cost method. Required for boost python*/
    DistanceType default_getDeletionCost(const TreeNodePtr node) const
    {  return NodeCost::getDeletionCost(node); }

    /** Returns the changing cost between /e i_node and /e r_node*/
    virtual DistanceType getChangingCost(const TreeNodePtr node1, const TreeNodePtr node2) const
    { 
	  if(bp::override f = this->get_override("getChangingCost"))
		return bp::call<DistanceType>(f.ptr(), treenode2python(node1), treenode2python(node2));
	  else return default_getChangingCost(node1,node2);
	}

    /** Default implementation of changing cost method. Required for boost python*/
    DistanceType default_getChangingCost(const TreeNodePtr node1, const TreeNodePtr node2) const
    {  return NodeCost::getChangingCost(node1, node2); }


    /** Returns the merging cost between /e i_node and /e r_node*/
  virtual DistanceType getMergingCost(const vector<TreeNodePtr> node1, const TreeNodePtr node2) const
    { 
	  if(bp::override f = this->get_override("getMergingCost"))
	    return bp::call<DistanceType>(f.ptr(), treenodelist2python(node1), treenode2python(node2)); //remplacer treenode2python par un vecteur de treenode ...
	  else return default_getMergingCost(node1,node2);
	}

    /** Default implementation of merging cost method. Required for boost python*/
  DistanceType default_getMergingCost(const vector<TreeNodePtr> node1, const TreeNodePtr node2) const
    {  return NodeCost::getMergingCost(node1, node2); }

    /** Returns the splitting cost between /e i_node and /e r_node*/
  virtual DistanceType getSplittingCost(const TreeNodePtr node1, const vector<TreeNodePtr> node2) const
    { 
	  if(bp::override f = this->get_override("getSplittingCost"))
	    return bp::call<DistanceType>(f.ptr(), treenode2python(node1), treenodelist2python(node2)); //remplacer treenode2python par un vecteur de treenode ...
	  else return default_getSplittingCost(node1,node2);
	}

    /** Default implementation of splitting cost method. Required for boost python*/
  DistanceType default_getSplittingCost(const TreeNodePtr node1, const vector<TreeNodePtr> node2) const
    {  return NodeCost::getSplittingCost(node1, node2); }

};

// A smart pointer on a PyNodeCost
typedef boost::shared_ptr<PyNodeCost> PyNodeCostPtr;

 

void export_NodeCost() {

	class_<PyNodeCost,PyNodeCostPtr,boost::noncopyable>
    ("NodeCost", init<>("Constructor of NodeCost")) 
    .def("getInsertionCost", &NodeCost::getInsertionCost, &PyNodeCost::default_getInsertionCost)
    .def("getDeletionCost", &NodeCost::getDeletionCost, &PyNodeCost::default_getDeletionCost)
    .def("getChangingCost", &NodeCost::getChangingCost, &PyNodeCost::default_getChangingCost)
    .def("getMergingCost", &NodeCost::getMergingCost, &PyNodeCost::default_getMergingCost)
    .def("getSplittingCost", &NodeCost::getSplittingCost, &PyNodeCost::default_getSplittingCost)
	;

}
