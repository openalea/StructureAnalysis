/*------------------------------------------------------------------------------
 *                                                                              
 *        VPlants.Stat_Tool : VPlants Statistics module
 *                                                                              
 *        Copyright 2006-2007 INRIA - CIRAD - INRA                      
 *                                                                              
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>         
 *                                                                              
 *        Distributed under the GPL 2.0 License.                               
 *        See accompanying file LICENSE.txt or copy at                          
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *                                                                              
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr                    
 *       
 *        $Id: scenecontainer_wrap.cpp 559 2007-05-25 12:25:30Z dufourko $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include <boost/python.hpp>

using namespace boost::python;


// Wrappers




// Overloads

// Format_error
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_update_overloads_1_3, update, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_correction_update_overloads_2_4, 
				       correction_update, 2, 4)



// Boost.Python Wrapper export function
void class_base()
{
  // _Format_error
  class_< Format_error >("_Format_error", init< optional< int > >())
    .def("init", &Format_error::init)
    .def("update", &Format_error::update, Format_error_update_overloads_1_3())
    .def("correction_update", (void (Format_error::*)(const char*, const char*, int, int) )
	 &Format_error::correction_update, Format_error_correction_update_overloads_2_4())
    .def("correction_update", (void (Format_error::*)(const char*, int, int, int) )
	 &Format_error::correction_update, Format_error_correction_update_overloads_2_4())
    .def("get_nb_error", &Format_error::get_nb_error)
    .def("get_max_nb_error", &Format_error::get_max_nb_error)
        .def(self_ns::str(self))
    ;

  

  // _Histogram
  class_<Histogram>("_Histogram", init<const Histogram&>());

}
