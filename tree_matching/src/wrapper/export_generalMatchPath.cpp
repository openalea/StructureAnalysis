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
//#include <hash_map>
#include "generalMatchPath.h"

using namespace boost::python;
#define bp boost::python



class PyGeneralMatchPath : public GeneralMatchPath, public bp::wrapper<GeneralMatchPath>
{
public:

    //Constructor
  PyGeneralMatchPath(boost::python::object in_list,  boost::python::object ref_list, boost::python::object in_successor, boost::python::object ref_predecessor) :  
      GeneralMatchPath( extract_vec<int>(in_list)(),
                        extract_vec<int>(ref_list)(),
                        extract_vec<std::vector<int>, extract_vec<int> >(in_successor)(),
                        extract_vec<std::vector<int>, extract_vec<int> >(ref_predecessor)()), 
      bp::wrapper<GeneralMatchPath>() 
  {   }

    //Constructor
  PyGeneralMatchPath(boost::python::object in_list,  boost::python::object ref_list, boost::python::object in_successor, boost::python::object ref_predecessor, boost::python::object in_capacities) :  
      GeneralMatchPath( extract_vec<int>(in_list)(),
                        extract_vec<int>(ref_list)(),
                        extract_vec<std::vector<int>, extract_vec<int> >(in_successor)(),
                        extract_vec<std::vector<int>, extract_vec<int> >(ref_predecessor)(),
                        extract_vec<int>(in_capacities)()), 
      bp::wrapper<GeneralMatchPath>() 
  {   }

  virtual ~PyGeneralMatchPath() {}


  /*  void py_make( boost::python::object in_list,  boost::python::object ref_list){
    NodeList input_list = extract_vec<int,boost::python::extract,std::vector<int> >(in_list)();
    NodeList reference_list = extract_vec<int,boost::python::extract,std::vector<int> >(ref_list)();
    this->make2(input_list,reference_list);
    }*/

    /** Returns the edge cost between /e i_node and /e r_node*/
  virtual DistanceType edgeCost(int node1, int node2) 
    { 
	  if(bp::override f = this->get_override("edgeCost")){
	    return bp::call<DistanceType>(f.ptr(), bp::object(node1), bp::object(node2));
      }
	  else return default_edgeCost(node1,node2);
	}

    /** Default implementation of edge cost method. Required for boost python*/
    DistanceType default_edgeCost(int node1, int node2) {
      return GeneralMatchPath::edgeCost(node1, node2); }


  boost::python::object py_minCostFlow(){
    VertexVector map_list;
    map_list.resize(this->nbVertex,-1);
    boost::python::list result;
    VertexVector::iterator it = map_list.begin();
    VertexVector::iterator begin = map_list.begin();
    VertexVector::iterator end = map_list.end();
    begin++;end--;
    int i = 1;
    //std::cerr<<"Mapping : "<<dist<<std::endl;
    DistanceType dist = this->minCostFlow(map_list);
    for (it = begin;it!=end;it++, i++){
      result.append(boost::python::make_tuple(this->who(i),who(*it)));
    }
    return boost::python::make_tuple(boost::python::object(dist),result);
  }
};

// A smart pointer on a PyGeneralMatchPath
typedef boost::shared_ptr<PyGeneralMatchPath> PyGeneralMatchPathPtr;
  

void export_GeneralMatchPath() {

	class_<PyGeneralMatchPath,PyGeneralMatchPathPtr,boost::noncopyable>
	  ("GeneralMatchPath", init<boost::python::object,
                                boost::python::object,
                                boost::python::object,
                                boost::python::object,
                                optional<boost::python::object> >("Constructor of GeneralMatchPath"))
	  .def( "bipartiteMatching", &PyGeneralMatchPath::py_minCostFlow,"Get the list of mapped vertices")
	  .def("edgeCost", &GeneralMatchPath::edgeCost, &PyGeneralMatchPath::default_edgeCost)
	;

}
