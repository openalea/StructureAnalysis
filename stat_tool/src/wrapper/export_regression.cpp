/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Dur&& <Jean-Baptiste.Dur&&@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_vectors.cpp 6274 2009-04-23 13:54:16Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/distribution.h"
#include "stat_tool/mv_mixture.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"
using namespace boost::python;
using namespace boost;




class RegressionKernelWrap
{
public:

};

void class_regression_kernel()
{
  class_< Regression_kernel  >("_Regression_kernel", "Regression kernel class")
    .def( init<int, int, int>())
    .def("get_ident", &Regression_kernel::get_ident,"return ident")
    .def("get_min_value", &Regression_kernel::get_min_value,"return min value")
    .def("get_max_value", &Regression_kernel::get_max_value,"return max value")
    .def("get_regression_df", &Regression_kernel::get_regression_df,"return regression df")
    .def("get_residual_df", &Regression_kernel::get_residual_df,"return residual df")
    .def("get_nb_parameter", &Regression_kernel::get_nb_parameter,"return nb parameter")
    .def("get_parameter", &Regression_kernel::get_parameter,python::args("index"),"return parameter")
    .def("get_point", &Regression_kernel::get_point,python::args("index"),"return point")

/*
Regression_kernel(const Regression_kernel &regression) { copy(regression); }
*/
  ;
}

class RegressionWrap
{
public:


  WRAP_METHOD_FILE_ASCII_WRITE(Regression);
  WRAP_METHOD_SPREADSHEET_WRITE(Regression);


  static double get_residual(Regression &input, int index)
  {
    double ret;
    ostringstream error_message;
    error_message << "index not in valid range" << endl;\
    CHECK(index, 0, input.get_nb_vector());
    ret = input.get_residual(index);
    return ret;
  }


};

void class_regression()
{
  class_< Regression, bases< STAT_interface > >
    ("_Regression", "Regression class")
    // Python Operators
    //
    .def(init <int, int, int, Vectors>())
    .def(init <Regression>())
    .def(self_ns::str(self)) // __str__
    .def("__len__", &Regression::get_nb_vector)  //__len__
    .def("get_nb_vector", &Regression::get_nb_vector, "Return nb_vector")
    .def("get_residual", RegressionWrap::get_residual, ARGS("int"),"Return nb_vector")
    .def("file_ascii_write", RegressionWrap::file_ascii_write, "Save regression summary into a file")
    .def("file_spreadsheet_write", RegressionWrap::spreadsheet_write, "Save regression summary into a CSV file")
    DEF_RETURN_VALUE_NO_ARGS("get_vectors", &Regression::get_vectors, "return vectors")
    ;

}




