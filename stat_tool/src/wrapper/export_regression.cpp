/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2014 CIRAD/INRA/Inria Virtual Plants
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
 *                        Jean-Baptiste Durand <Jean-Baptiste.Durand@imag.fr>
 *                        Samuel Dufour-Kowalski <samuel.dufour@sophia.inria.fr>
 *                        Christophe Pradal <christophe.prada@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_regression.cpp 6274 2009-04-23 13:54:16Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/



#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/multivariate_mixture.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;


#define WRAP RegressionKernelWrap

class RegressionKernelWrap
{

public:

  static boost::shared_ptr<RegressionKernel>
  regression_constructor(int ident, int min_value, int max_value)
  {
    RegressionKernel *ret = NULL;
    ret = new RegressionKernel(ident, min_value, max_value);

    return boost::shared_ptr<RegressionKernel>(ret);
  }

};



void class_regression_kernel()
{
  class_< RegressionKernel>("_RegressionKernel", "Regression kernel class")
    //.def( init<int, int, int>()) // seems to cause memory problem inside python
    .def("__init__", make_constructor(WRAP::regression_constructor),  "Build from 3 ints" )

/*    .add_property("ident", &RegressionKernel::get_ident, "returns ident")
    .add_property("min_value", &RegressionKernel::get_min_value, "returns min value")
    .add_property("max_value", &RegressionKernel::get_max_value, "returns max value")
    .add_property("regression_df", &RegressionKernel::get_regression_df, "returns regression df")
    .add_property("residual_df", &RegressionKernel::get_residual_df, "returns residual df")
    .add_property("nb_parameter", &RegressionKernel::get_nb_parameter, "returns nb parameter")

    .def("get_parameter", &RegressionKernel::get_parameter,args("index"), "returns parameter")
    .def("get_point", &RegressionKernel::get_point,args("index"), "returns point") */

    //DONE
    //TODO: check index validity of get_parameter and get_point
    /*
    RegressionKernel(const RegressionKernel &regression) { copy(regression); }
    */
  ;
}
#undef WRAP




#define WRAP RegressionWrap
class WRAP
{
public:

  WRAP_METHOD_FILE_ASCII_WRITE( Regression);
  WRAP_METHOD_SPREADSHEET_WRITE( Regression);

  static double
  get_residual(Regression &input, int index)
  {
    double ret;
    ostringstream error_message;
    error_message << "index not in valid range" << endl;
    CHECK(index, 0, input.get_nb_vector());
    ret = input.get_residual(index);
    return ret;
  }

};

void class_regression()
{
  class_< Regression, bases<StatInterface> >
  ("_Regression", "Regression class")
    .def(init <int, int, int, const Vectors&>())
    .def(init <const Regression &>())

    .def(self_ns::str(self)) // __str__
    .def("__len__", &Regression::get_nb_vector)  //__len__

    .add_property("nb_vector", &Regression::get_nb_vector, "Return nb_vector")

    .def("get_residual", WRAP::get_residual, args("int"),"Return nb_vector")
    .def("file_ascii_write", WRAP::file_ascii_write, "Save regression summary into a file")
    .def("file_spreadsheet_write", WRAP::spreadsheet_write, "Save regression summary into a CSV file")

    DEF_RETURN_VALUE_NO_ARGS("get_vectors", &Regression::get_vectors, "return vectors")
    ;

  //DONE
  /*
    bool plot_write(StatError &error , const char *prefix , const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;
    */
}
#undef WRAP




