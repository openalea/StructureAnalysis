/*------------------------------------------------------------------------------
 *
 *        VPlants.Stat_Tool : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Thomas.Cokelaer@inria.fr
 *
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


#include "stat_tool/stat_tools.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"


#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;


void class_non_parametric_process() {

        class_<Nonparametric_process >
        ("_Nonparametric_process", "Non parametric process")

        .def(init< optional< int, int, int> > ())
        .def("get_nb_value", &Nonparametric_process::get_nb_value, "returns number of values")
        .def("get_nb_state", &Nonparametric_process::get_nb_state, "returns number of states")

       ;
        /*
   Nonparametric_process(int inb_state , int inb_value , double **observation_probability);
   Nonparametric_process(int inb_state , Distribution **pobservation);
   Nonparametric_process(const Nonparametric_process &process , char manip = 'c' , int state = I_DEFAULT);

   int nb_parameter_computation(double min_probability) const;


   Distribution* get_observation(int state) const   { return observation[state]; }

   std::ostream& ascii_print(std::ostream &os, Histogram **empirical_observation,
                 bool exhaustive, bool file_flag) const;

   std::ostream& spreadsheet_print(std::ostream &os, Histogram **empirical_observation=NULL) const;

   bool plot_print(const char *prefix, const char *title, int process, Histogram **empirical_observation = NULL) const;


*/
}
