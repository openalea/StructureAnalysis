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

#include "treenode.h"
#include "extract_list.h"
#include "export_list.h"
#include <boost/python.hpp>
using namespace boost::python;
#define bp boost::python

void py_setchildlist(TreeNode* tnode, object o){
  tnode -> setChildList(extract_vec<int>(o)());
}

object py_getchildlist(TreeNode* tnode){
  return make_list(tnode->getChildList())();
}

std::string treenode_str( TreeNode* tnode ) 
{ 
  stringstream ss; 
  ss<<"ID \t : "<<(int)tnode->getId()<<endl;
  ss<<"FATHER \t : "<<(int)tnode->father()<<endl;
  ss<<"DEPTH \t : "<<(int)tnode->depth()<<endl;
  ss<<"CHILD LIST \t : [ ";
  for (int i=0;i<tnode->getChildNumber()-1;i++)
    ss<<tnode->getChild(i)<<" , ";
  if (tnode->getChildNumber()>0)
    ss<<tnode->getChild(tnode->getChildNumber()-1)<<" ] "<<endl;
  else
    ss<<"]"<<endl;
  for (int i=0;i<tnode->getValueSize();i++)
    ss<<"ARG["<<i<<"] \t : "<<tnode->getValue(i)<<endl;
  return ss.str(); 
} 

std::string treenode_repr( TreeNode* tnode ) 
{ 
  stringstream ss; 
  ss<<"<TreeNode(id="<< tnode->getId();
  if( tnode->father() != -1)
	ss<<",father_id="<<tnode->father();
  ss<<",depth="<<(int)tnode->depth();
  ss<<",children=[";
  for (int i=0;i<tnode->getChildNumber();++i){
    if (i>0) ss << ",";
    ss<<tnode->getChild(i)<<",";
  }
  ss<<"],values=[";
  for (int i=0;i<tnode->getValueSize();++i){
    if (i>0) ss << ",";
    ss <<tnode->getValue(i);
  }
  ss<<"]) at 0x" << tnode << ">";
  return ss.str(); 
} 

void export_TreeNode() {

  class_<TreeNode,TreeNodePtr,boost::noncopyable>
	  ("TreeNode", init<int,int>("TreeNode(id,father_id)",(bp::arg("id"),bp::arg("father_id")=-1)))
    .add_property("father",&TreeNode::father,&TreeNode::putFather)
    .add_property("id",&TreeNode::getId,&TreeNode::putId)
    .def( "putValue", &TreeNode::putValue,"Set the Values",(bp::arg("index"),bp::arg("new_value")=0.))
    .def( "addValue", &TreeNode::addValue,"Set the Values",(bp::arg("new_value")=0.))
    .def( "addChild", &TreeNode::addChild,"Put a child id in the Child list",(bp::arg("child_id")))
    .def( "getChild", &TreeNode::getChild,"Get a child id from his position id in the Child list")
    .def( "setChildList", &py_setchildlist,"Set the ChildList",(bp::arg("child_list")))
    .def( "getChildList", &py_getchildlist,"Get the ChildList")
    .def( "getChildListSize", &TreeNode::getChildNumber,"Id of Children")
    .def( "__repr__", treenode_repr )
    .def( "__str__", treenode_str )
	;

}
