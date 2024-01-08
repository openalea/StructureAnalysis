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
#include "matching.h"
#include "treegraph.h"
using namespace boost::python;
#define bp boost::python

object py_getDistanceTable(Matching* m){
  return make_list<DistanceVectorTable,list_converter<DistanceVectorTable::value_type> >(m->getDistanceTable())();
}

boost::python::object py_getList(Matching * m, int i_tree,int r_tree)
{
  Sequence list_vtx;  
  m->getList(i_tree,r_tree,&list_vtx);
  list_vtx.reset();
  boost::python::list result;
  do {
    result.append(boost::python::make_tuple(list_vtx.getCurrent()->getIV(),
			     list_vtx.getCurrent()->getRV(),
			     list_vtx.getCurrent()->getCost()));
  } while(list_vtx.next());
  return result;
}



void export_Matching() {

  class_<Matching>
    ("Matching", init<TreeGraphPtr, TreeGraphPtr, NodeCostPtr, int>("Matching(TreeGraph tree1, TreeGraph tree2, NodeCost nodecost, MDTableType tabletype = 0)",
     (bp::arg("tree1"),bp::arg("tree2"),bp::arg("nodecost"),bp::arg("tabletype")=STD)))
    .def( "match", &Matching::match,"Comparison of tree graph")
    .def( "__call__", &Matching::match,"Comparison of tree graph")
    .def( "getDBT", &Matching::getDBT,"Get Distance Between Trees",(bp::arg("index1"))=0,(bp::arg("index2")=0))
    .def( "getDistanceTable", &py_getDistanceTable,"Get DistanceTable Between Trees")
    .def( "getList", &py_getList,"Get Matching  list between Trees",(bp::arg("index1"))=0,(bp::arg("index2")=0))
    .def( "getChoiceTable", &Matching::getChoiceTable,"Get Choice List",return_internal_reference<>())
    .add_property("verbose",make_getter(&Matching::verbose),make_setter(&Matching::verbose))
	;

}

void export_ExtMatching() {

  class_<ExtMatching,bases<Matching> >
    ("ExtMatching", init<TreeGraphPtr, TreeGraphPtr, NodeCostPtr, int>("ExtMatching(TreeGraph, TreeGraph, NodeCost, MDTableType)"));
}
