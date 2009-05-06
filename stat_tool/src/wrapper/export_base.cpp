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
 *        $Id$
 *
 *-----------------------------------------------------------------------------*/

#include "export_base.h"
#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"
using namespace boost::python;
using namespace boost;




// Overloads

// Format_error
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_update_overloads_1_3, update, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_correction_update_overloads_2_4,
				       correction_update, 2, 4)


// Boost.Python Wrapper export function

void class_constant()
{
  //constant
  scope().attr("I_DEFAULT") = I_DEFAULT;
  scope().attr("D_DEFAULT") = D_DEFAULT;
  scope().attr("D_INF") = D_INF;
  scope().attr("MAX_DIFF_BOUND") = MAX_DIFF_BOUND;
  scope().attr("MAX_MEAN") = MAX_MEAN;

  enum_<stat_tool::wrap_util::UniqueInt<6, 0> >("VariableTypeBis")
    .value("INT_VALUE", INT_VALUE)
    .value("REAL_VALUE", REAL_VALUE)
    .value("STATE", STATE)
    .value("OLD_INT_VALUE", OLD_INT_VALUE)
    .value("NB_INTERNODE", NB_INTERNODE)
    .value("AUXILIARY", AUXILIARY)
    .export_values()
  ;


}

// Format Error

void class_format_error()
{
  // _Format_error
  class_< Format_error >("_FormatError", init< optional< int > >())
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

}

class StatInterfaceWrap
{
public:

  WRAP_METHOD_ASCII_WRITE(STAT_interface);
  WRAP_METHOD_PLOT_WRITE(STAT_interface);
  WRAP_METHOD_SPREADSHEET_WRITE(STAT_interface);
};



void class_stat_interface()
{
  class_< STAT_interface, boost::noncopyable > ("_StatInterface", no_init)
    .def("ascii_write", &StatInterfaceWrap::ascii_write, args("exhaustive"),
    	"Return a string containing the object description (exhaustive or not)")

    .def("plot_write", &StatInterfaceWrap::plot_write,	args("prefix", "title"),
  		"Write GNUPLOT files (with prefix)")

    .def("spreadsheet_write", &StatInterfaceWrap::spreadsheet_write, args("filename"),
   		"Write object to filename (spreadsheet format)")

    .def("get_plotable", &STAT_interface::get_plotable, return_value_policy< manage_new_object >(),
  		"Return a plotable object" )
    ;

}


void class_forward()
{
  class_<Forward , boost::noncopyable ,bases<Parametric> > ("_Forward", no_init)

;
/*
    Forward(int inb_value = 0 , int iident = NONPARAMETRIC ,  int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT , double iparameter = D_DEFAULT , double iprobability = D_DEFAULT)    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Forward(const Parametric &dist , int ialloc_nb_value = I_DEFAULT) :Parametric(dist , 'c' , ialloc_nb_value) { computation(dist); }
    Forward(const Forward &forward , int ialloc_nb_value = I_DEFAULT)  :Parametric((Parametric&)forward , 'c' , ialloc_nb_value) {}

    void computation(const Parametric &dist);
*/
}


