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
    .def("mtg_write",(bool (TreeGraph::*)( char *))&TreeGraph::mtg_write,"Save TreeGraph as MTG file",(bp::arg("path")))
    .def( "__repr__", treegraph_str )
    .def( "__str__", treegraph_str )
    .def( "__repr__", treegraph_str )
	;

}
