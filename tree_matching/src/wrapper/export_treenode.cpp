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


boost::python::object  py_getValue(TreeNode * n, size_t index)
{
	boost::any o = n->getAnyValue(index);
	if (o.type() == typeid(double)) return boost::python::object(boost::any_cast<double>(o));
	else if (o.type() == typeid(int)) return boost::python::object(boost::any_cast<int>(o));
	else return boost::any_cast<boost::python::object>(o);
}

void py_setValue(TreeNode * n, size_t index, boost::python::object val)
{ n->setTypedValue(index,val); }

void  py_appendValue(TreeNode * n, boost::python::object val)
{ n->appendTypedValue<boost::python::object>(val); }

boost::python::object  py_getValues(TreeNode * n)
{ 
	boost::python::list l;
	for(TreeNode::ValueVector::const_iterator it = n->getValueList().begin(); 
		it != n->getValueList().end(); ++it)
		l.append(boost::any_cast<boost::python::object>(*it));
	return l; 
}

std::string treenode_str( TreeNode* tnode ) 
{ 
  stringstream ss; 
  ss<<"ID \t : "<<(int)tnode->getId()<<endl;
  ss<<"FATHER \t : "<<(int)tnode->getFather()<<endl;
  ss<<"DEPTH \t : "<<(int)tnode->getDepth()<<endl;
  ss<<"CHILD LIST \t : [ ";
  if(tnode->hasChild()){
	  for (size_t i=0; i < tnode->getChildNumber();i++){
		if (i>0) ss <<" , ";
		ss << tnode->getChild(i);
	  }
  }
  ss<<"]"<<endl;
  for (int i=0;i<tnode->getValueSize();i++) {
    ss << "ARG["<<i<<"] \t : ";
	ss << extract<std::string>(str(py_getValue(tnode,i)))() ;
	ss << endl;
  }
  return ss.str(); 
} 

std::string treenode_repr( TreeNode* tnode ) 
{ 
  stringstream ss; 
  ss<<"<TreeNode(id="<< tnode->getId();
  if( tnode->getFather() != -1)
	ss<<",father_id="<<tnode->getFather();
  ss<<",depth="<<(int)tnode->getDepth();
  ss<<",children=[";
  for (int i=0;i<tnode->getChildNumber();++i){
    if (i>0) ss << ",";
    ss<<tnode->getChild(i)<<",";
  }
  ss<<"],values=[";
  for (int i=0;i<tnode->getValueSize();++i){
    if (i>0) ss << ",";
	ss << extract<std::string>(str(py_getValue(tnode,i)))();
  }
  ss<<"]) at 0x" << tnode << ">";
  return ss.str(); 
} 

static boost::python::object * TreeNodeBuilder = NULL;

TreeNodePtr py_factory_build(int id, int father) {
  if (TreeNodeBuilder == NULL) return TreeNodePtr();
	return call<TreeNodePtr>((*TreeNodeBuilder).ptr(),id,father);
}

void py_factory_setBuilder(TreeNode::Factory * factory, boost::python::object b){
	TreeNodeBuilder = new boost::python::object(b);
	factory->setBuilder(&py_factory_build);
}

void py_factory_finalize() {
  if (TreeNodeBuilder == NULL) delete TreeNodeBuilder;
}

void export_TreeNode() {

  scope tn = class_<TreeNode,TreeNodePtr,boost::noncopyable>
	  ("TreeNode", init<int,int>("TreeNode(id,father_id)",(bp::arg("id"),bp::arg("father_id")=-1)))
    .add_property("father",&TreeNode::getFather,&TreeNode::setFather)
    .add_property("id",&TreeNode::getId,&TreeNode::setId)

    .def( "putValue", &TreeNode::getTypedValue<DistanceType>,"Set the Values",(bp::arg("index"),bp::arg("new_value")=0.))
    .def( "addValue", &TreeNode::appendTypedValue<DistanceType>,"Set the Values",(bp::arg("new_value")=0.))

    .def( "__getitem__", &py_getValue,"Get a value",(bp::arg("new_value")=0.))
    .def( "__setitem__", &py_setValue,"Set a value",(bp::arg("new_value")=0.))
    .def( "append", &py_appendValue,"Append a value",(bp::arg("new_value")=0.))
	.add_property("values",&py_getValues)

    .def( "addChild", &TreeNode::addChild,"Put a child id in the Child list",(bp::arg("child_id")))
    .def( "getChild", &TreeNode::getChild,"Get a child id from his position id in the Child list")

    .def( "setChildList", &py_setchildlist,"Set the ChildList",(bp::arg("child_list")))
    .def( "getChildList", &py_getchildlist,"Get the ChildList")
    .def( "getChildListSize", &TreeNode::getChildNumber,"Number of Children")

    .def( "factory", &TreeNode::factory,return_value_policy<reference_existing_object>())
	.staticmethod("factory")

    .def( "__repr__", treenode_repr )
    .def( "__str__", treenode_str )
	;

  class_<TreeNode::Factory,boost::noncopyable>("Factory",no_init)
	  .def("setBuilder",&py_factory_setBuilder,args("builder"))
	  .def("setBuilderToDefault",&TreeNode::Factory::setBuilderToDefault)
	  .def("__call__",&TreeNode::Factory::build)
	  .def("build",&TreeNode::Factory::build)
	  ;
}
