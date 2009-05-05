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
 *        $Id: export_tops.cpp 6169 2009-04-01 16:42:59Z cokelaer $
 *
 *-----------------------------------------------------------------------------*/
#include "wrapper_util.h"

#include <math.h>
#include <iomanip>
#include <sstream>
#include <tool/rw_tokenizer.h>
#include <tool/rw_cstring.h>
#include <tool/rw_locale.h>


#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/sequence_label.h"
#include "tool/config.h"


#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
//using namespace stat_tool;

////////////////////// Export tops ////////////////////////////////////////

class FunctionWrap {

public:

  static boost::shared_ptr<Function> function_from_list(int ident, int length, boost::python::list& input_parameter)
  {
    Function *ret = NULL;
    double *parameter;

    int nb_param = len(input_parameter);
    parameter = new double[nb_param];

    if (nb_param == 0)
    {
        sequence_analysis::wrap_util::throw_error("input list cannot be empty");
    }

    for (int i=0; i<nb_param; i++)
    {
        parameter[i] = boost::python::extract<double> (input_parameter[i]);
    }

    ret = new Function(ident, length, parameter);
    if (!ret)
        sequence_analysis::wrap_util::throw_error("error while calling Function constructor.");

    return boost::shared_ptr<Function>(ret);
  }

};

// Boost declaration

void class_function() {

  class_<Function, bases<Regression_kernel> > ("_Function", "Function")
    .def("__init__", make_constructor(FunctionWrap::function_from_list))
    .def(init<int, int>())
    .def(init<const Function &>())
    .def("get_residual",&Function::get_residual,args("index"), "get residual")
    .def("get_frequency",&Function::get_frequency,args("index"), "get frequency")

    ;
  //DONE
  //TODO : check validity of the input index on get_residual and get_frequency
}
