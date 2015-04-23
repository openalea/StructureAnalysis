/*------------------------------------------------------------------------------
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2014 CIRAD/INRA/Inria Virtual Plants
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id:  $
 *
 *----------------------------------------------------------------------------*/



#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
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
using namespace stat_tool;
using namespace sequence_analysis;



#define WRAP CorrelationWrap
class CorrelationWrap {

public:

  // Merge
  static Correlation*
  merge(const Correlation& input, const boost::python::list& input_cor)
  {
    CREATE_ARRAY(input_cor, const Correlation *, data);
    SIMPLE_METHOD_TEMPLATE_1(input, merge, Correlation,
        data_size, data.get());
  }


  static bool
  white_noise_correlation_order(Correlation& input, int order)

  {
    StatError error;
    bool ret;
    ret = input.white_noise_correlation(error, order);
    FOOTER;
  }

  static bool
  white_noise_correlation_filter(Correlation& input,
		  boost::python::list input_filter)

  {
   StatError error;
   bool ret;
   int nb_point = len(input_filter);
   double sum=0;
   int i=0;
 
   sequence_analysis::wrap_util::auto_ptr_array<double> filter(new double[nb_point * 2 + 1]);
    

   nb_point--;
   for (i = 0; i < nb_point; i++)
   {
     filter[i] = boost::python::extract<double> (input_filter[i]);
     filter[2*nb_point -i] = filter[i];
     sum += 2 * filter[i];
   }
   //i = n
   filter[i] = extract<double> (input_filter[i]);
   sum += filter[i];

   for (i = 0; i < 2 * nb_point + 1; i++)
   {
     filter[i] = filter[i] / sum;
   }

    
   ret = input.white_noise_correlation(error, 2 * nb_point +1, filter.get());
   //if (!ret) throw error
  
  }

  static bool
  white_noise_correlation_dist(Correlation& input,
 		  const Distribution& dist)
  {
    StatError error;
    bool ret;
    ret = input.white_noise_correlation(error, dist);
    //if (!ret) throw error
   }


};


void class_correlation() {

  class_<Correlation, bases<StatInterface> >
  ("_Correlation", "Correlation")
    .def(init<int, int, int, int>())
    .def(init<int, int, bool, int>())
    .def(init<Correlation>())

    .def(self_ns::str(self)) //__str__

    .add_property("type", &Correlation::get_type)

    .def("get_variable_type", &Correlation::get_variable_type, args("index"))
    .def("get_variable1", &Correlation::get_variable1, args("index"))
    .def("get_variable2", &Correlation::get_variable2, args("index"))
    .def("get_white_noise", &Correlation::get_white_noise, args("lag"))

    .def("white_noise_correlation_dist", WRAP::white_noise_correlation_dist, args("dist"), "todo")
    .def("white_noise_correlation_order", WRAP::white_noise_correlation_order, args("order"), "todo")
    .def("white_noise_correlation_filter", WRAP::white_noise_correlation_filter, args("filter"), "todo")

    DEF_RETURN_VALUE("merge", WRAP::merge, args("list"),"todo")
    ;
}



#undef WRAP
