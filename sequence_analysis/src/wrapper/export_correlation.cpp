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


#include "stat_tool/stat_tools.h"
#include "stat_tool/vectors.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/sequence_label.h"
#include "tool/config.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace sequence_analysis;

#define WRAP CorrelationWrap
class CorrelationWrap {

public:

  // Merge
  static Correlation*
  merge(const Correlation& input_cor, const boost::python::list& correlations)
  {
    int nb_cor = len(correlations);
    sequence_analysis::wrap_util::auto_ptr_array<const Correlation *> sequens(
        new const Correlation*[nb_cor]);
    for (int i = 0; i < nb_cor; i++)
      sequens[i] = extract<Correlation*> (correlations[i]);

    SIMPLE_METHOD_TEMPLATE_1(input_cor, merge, Correlation,
        nb_cor, sequens.get());
  }


  static bool
  white_noise_correlation_order(Correlation& input, int order)

  {
	  Format_error error;
	  bool ret;
	  ret = input.white_noise_correlation(error, order);
	  //if (!ret) throw error
  }

  static bool
  white_noise_correlation_filter(Correlation& input,
		  boost::python::list input_filter)

  {
   Format_error error;
   bool ret;
   int nb_point = len(input_filter);
   double *filter;
   filter = new double[2*nb_point+1];

   for (int i = 0; i < nb_point; i++)
   {
     filter[i] = boost::python::extract<double> (input_filter[i]);
     filter[2*nb_point -i] = filter[i];
   }
   filter[nb_point] = 0; //todo check that!!

   ret = input.white_noise_correlation(error, 2 * nb_point +1, filter);
   //if (!ret) throw error
  }

  static bool
  white_noise_correlation_dist(Correlation& input,
 		  const Distribution& dist)
  {
    Format_error error;
    bool ret;
    ret = input.white_noise_correlation(error, dist);
    //if (!ret) throw error
   }


};


void class_correlation() {

  class_<Correlation, bases<STAT_interface> >
  ("_Correlation", "Correlation")
    .def(init<int, int, int, int>())
    .def(init<int, int, bool, int>())
    .def(init<Correlation>())

    .def(self_ns::str(self)) //__str__

    .add_property("get_type", &Correlation::get_type)

    .def("get_variable_type", &Correlation::get_variable_type,args("index"))
    .def("get_variable1", &Correlation::get_variable1,args("index"))
    .def("get_variable2", &Correlation::get_variable2,args("index"))
    .def("get_white_noise", &Correlation::get_white_noise,args("lag"))
    DEF_RETURN_VALUE("merge", WRAP::merge, args("list"),"todo")
    .def("white_noise_correlation_dist", WRAP::white_noise_correlation_dist, args("dist"), "todo")
    .def("white_noise_correlation_order", WRAP::white_noise_correlation_order, args("order"), "todo")
    .def("white_noise_correlation_filter", WRAP::white_noise_correlation_filter, args("filter"), "todo")


    ;

//todo
/*
  std::ostream& line_write(std::ostream &os) const;
*/
}



#undef WRAP
