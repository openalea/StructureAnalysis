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


#include "extract_list.h"
#include "export_list.h"
#include <boost/python.hpp>
#include "matching.h"
#include "treegraph.h"
using namespace boost::python;
#define bp boost::python

object py_getDistanceTable(Matching* m){
  return make_list<DistanceTable,list_converter<DistanceTable::value_type> >(m->getDistanceTable())(); //ou value_type ???
}

 

void export_Matching() {

  class_<Matching>
    ("Matching", init<TreeGraph, TreeGraph, NodeCost>("Matching(TreeGraph, TreeGraph, NodeCost)"))
    .def( "match", &Matching::match,"Comparison of tree graph")
    .def( "getDBT", &Matching::getDBT,"Get Distance Between Trees",(bp::arg("index")),(bp::arg("index")))
    .def( "getDistanceTable", &py_getDistanceTable,"Get DistanceTable Between Trees")
	;

}
