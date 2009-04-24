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
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/tops.h"


#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;




////////////////////// Export tops ////////////////////////////////////////

class TopParametersWrap
{

public:

  static boost::shared_ptr<Top_parameters> top_parameters_from_file(char* filename, int max_position)
  {
    Format_error error;
    Top_parameters *top_parameters = NULL;

    top_parameters = top_parameters_ascii_read(error, filename, max_position);

/*    if(!top_parameters)
    {
	  stat_tool::wrap_util::throw_error(error);
    }
*/
    return boost::shared_ptr<Top_parameters>(top_parameters);
  }

};



// Boost declaration

void class_top_parameters()
{

  class_< Top_parameters, bases< STAT_interface > >
    ("_Top_parameters", "Top parameters")
    .def("__init__", make_constructor(TopParametersWrap::top_parameters_from_file))
    .def(self_ns::str(self)) //__str__
    .def_readonly("get_probability", &Top_parameters::get_probability)
    .def_readonly("get_axillary_probability", &Top_parameters::get_axillary_probability)
    .def_readonly("get_rhythm_ratio", &Top_parameters::get_rhythm_ratio)
    .def("get_max_position", &Top_parameters::get_max_position)


    ;

}



class TopsWrap
{

public:

  static boost::shared_ptr<Tops> tops_from_file(char* filename, bool old_format)
  {
    Format_error error;
    Tops *tops = NULL;

    tops = tops_ascii_read(error, filename, old_format);

    if(!tops)
    {
	  stat_tool::wrap_util::throw_error(error);
    }

    return boost::shared_ptr<Tops>(tops);
  }


//     stat_tool::wrap_util::auto_ptr_array<const Parametric *>
//       dist(new const Parametric*[nb_dist]);

//     for(int i=0; i<nb_dist; i++)
//     dist[i] = boost::python::extract< Parametric *>(dists[i]);

//     conv = convolution_building(error, nb_dist, dist.get());



  static Histogram* get_nb_internode(const Tops& top)
  {

	  Histogram *h = NULL;

	  h = top.get_nb_internode();

	  return h;

  }


};


// Boost declaration

void class_tops()
{

  class_< Tops, bases< Sequences > >
    ("_Tops", "Tops")
    .def("__init__", make_constructor(TopsWrap::tops_from_file))
    .def(init<Sequences&>())


    .def(self_ns::str(self)) //__str__

    .def("get_max_position", &Tops::get_max_position)
    .def("get_nb_internode", TopsWrap::get_nb_internode,
    		return_value_policy<manage_new_object> ())
    .def("extract", &Tops::extract,
    		return_value_policy<manage_new_object> ())
    .def("shift", &Tops::shift,
    		return_value_policy<manage_new_object> ())

    ;
}
