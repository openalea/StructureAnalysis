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

#include "choicetable.h"
#include "extract_list.h"
#include "export_list.h"
#include <boost/python.hpp>
#include "treenode.h"
using namespace boost::python;
#define bp boost::python


boost::python::object py_getList(ChoiceTable * m, int i_tree,int r_tree, TreeGraphPtr T1, TreeGraphPtr T2)
{
  MatchRecordList list_vtx;  
  m->getList(i_tree,r_tree,T1, T2,&list_vtx);
  boost::python::list result;
  for(MatchRecordList::const_iterator it = list_vtx.begin(); it != list_vtx.end(); ++it)
	  result.append(make_tuple(it->first,it->second));
  return result;
}


void export_ChoiceTable() {

	class_<ChoiceTable> ("ChoiceTable", no_init)
    .def( "getList", &py_getList,"Get Matching  list between Trees",(bp::arg("index1"),bp::arg("index2"),bp::arg("tree1"),bp::arg("tree2")))
    .def( "dump", &ChoiceTable::dump, "Dump the choice table in a file",(bp::arg("filename"))) 
    .def( "load", &ChoiceTable::load, "Load a ChoiceTable from a file" ,(bp::arg("filename")))
	.staticmethod("load")
	;

}
