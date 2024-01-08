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

#include "treegraph.h"
#include "extract_list.h"
#include "export_list.h"
#include <boost/python.hpp>
#include "treenode.h"
using namespace boost::python;
#define bp boost::python


std::string treegraph_str( TreeGraph* tree ) 
{ 
  stringstream ss;
  ss<<"TREE SIZE\t: "<<(int)tree->getNbVertex()<<endl;
  ss<<"DEGREE\t = "<<(int)tree->getDegree()<<endl;
  for (int i = 0; i < tree->getNbVertex() ; i++)
    tree->getNode(i)->print();

//   TreeNode* tnode = tree->getTreeNode(tree->getRoot());
//   ss<<treenode_str(tnode);
  return ss.str(); 
} 


int py_next_iter(TreeGraph::const_pre_order_iterator * iter) {
    if (iter->atEnd()) {
         PyErr_SetString(PyExc_StopIteration, "index out of range");
         boost::python::throw_error_already_set();
    }
    int value = *(*iter);
    ++(*iter);
    return value;
}

boost::python::object py_subtree_list(TreeGraph * tg, int node){
    boost::python::list result;
    for(TreeGraph::const_pre_order_iterator it = tg->subtree_iterator_begin(node); it != tg->subtree_iterator_end(node); ++it)
        result.append(*it);
    return result;
}


inline void nullfunc(TreeGraph::const_pre_order_iterator * ) { }

void export_TreeGraph() {

	class_<TreeGraph,TreeGraphPtr,boost::noncopyable>
    ("TreeGraph", init<>("TreeGraph()"))
    .def( "addNode", (void (TreeGraph::*)(int,int))&TreeGraph::addNode,"Add a new node",(bp::arg("id"),bp::arg("father")=-1))
    .def( "addNode", (void (TreeGraph::*)(TreeNodePtr))&TreeGraph::addNode,"Add a new node",bp::args("node"))
    .def( "addValue", &TreeGraph::addValue,"Add a new value to an existing node",(bp::arg("index")),(bp::arg("new_value")))
    .def( "getNode", &TreeGraph::getNode,"Return the node referenced by id",(bp::arg("vertex")))
    .def( "father", &TreeGraph::father,"Return the father of vertex",(bp::arg("vertex")))
    .def( "child", &TreeGraph::child,"Return the id of child_number th child of node",(bp::arg("node")),(bp::arg("child_number")))
    .def( "getNbChild", &TreeGraph::getNbChild,"Get the number of children",(bp::arg("vertex")))
    .def( "getNbVertex", &TreeGraph::getNbVertex, "Get the number of vertex of the tree" )
    .def( "getNbDescendant", &TreeGraph::getNbDesc, "Get the number of descendant of a vertex",(bp::arg("vertex")) )
    .def( "getRoot", &TreeGraph::getRoot, "Get the root of the tree graph" )
    .def( "getDegree", &TreeGraph::getDegree, "Get the degree of the tree graph" )
    .def( "getDepth", &TreeGraph::getDepth, "Get the depth of the tree graph" )
    .def( "isLeaf", &TreeGraph::isLeaf, "Return is a vertex is a leaf",(bp::arg("vertex")) )
    .def( "empty", &TreeGraph::isNull, "Return is the tree graph is empty or not." )
    .def("mtg_write",(bool (TreeGraph::*)( char *))&TreeGraph::mtg_write,"Save TreeGraph as MTG file",(bp::arg("path")))
    .def( "__repr__", treegraph_str )
    .def( "__str__", treegraph_str )
    .def( "__repr__", treegraph_str )
    .def( "subtree_iter", &TreeGraph::subtree_iterator_begin)
    .def( "subtree_list", &py_subtree_list)
	;

	class_<TreeGraph::const_pre_order_iterator>("TreeGraphPreOrderIterator",no_init)
    .def("next",&py_next_iter)
    .def("__next__",&py_next_iter)
    .def("atEnd",&TreeGraph::const_pre_order_iterator::atEnd)
    .def("__iter__",&nullfunc, return_self<>())
    ;
}
