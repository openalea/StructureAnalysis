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

#include "wrapper_util.h"
#include "export_distribution.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"


#include <boost/python.hpp>
// definition of boost::python::len
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::make_constructor
#include <boost/python/make_constructor.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;
using namespace boost;


////////////////////// Export Parametric_model ////////////////////////////////




// Wrapper class

class ParametricModelWrap
{

public:

  static boost::shared_ptr<Parametric_model> parametric_model_from_file(char* filename)
  {
    Format_error error;
    Parametric_model *model = NULL;
    model = parametric_ascii_read(error, filename);

    if(!model) stat_tool::wrap_util::throw_error(error);
    return boost::shared_ptr<Parametric_model>(model);
  }


  static Distribution_data* simulation(const Parametric_model& p, int nb_element)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = p.simulation(error, nb_element);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }
  
  static Distribution_data* extract_data(const Parametric_model& p)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = p.extract_data(error);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static std::string survival_ascii_write(const Parametric_model& p)
  {
    std::stringstream s;
    std::string res;
    
    p.survival_ascii_write(s);
    res = s.str();
    
    return res;
  }


  
  static void survival_spreadsheet_write(const Parametric_model& p,
						const std::string& filename)
  {
    Format_error error;

    if(!p.survival_spreadsheet_write(error, filename.c_str()))
      stat_tool::wrap_util::throw_error(error);
      
  }


   
  static void survival_plot_write(const Parametric_model& p,
				  const std::string& prefix, const std::string& title)
  {
    Format_error error;

    if(!p.survival_plot_write(error, prefix.c_str(), title.c_str()))
      stat_tool::wrap_util::throw_error(error);
  }

  
  static MultiPlotSet* survival_get_plotable(const Parametric_model& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }

  
//   static void plot_write(const Parametric_model& p,
// 			 const std::string& prefix, const std::string& title,
// 			 const boost::python::list& dist_list)
//   {
//     Format_error error;

//     int nb_dist = boost::python::len(dist_list);
//     stat_tool::wrap_util::auto_ptr_array<const Distribution *> 
//       dists(new const Distribution*[nb_dist]);

//     const Distribution &d = (const Distribution&)(p);

//     if(!d.plot_write(error, prefix.c_str(), nb_dist, dists.get(), title.c_str()))
//       stat_tool::wrap_util::throw_error(error);
//   }


  static MultiPlotSet* get_plotable(const Parametric_model& p, 
				const boost::python::list& dist_list)
  {
    Format_error error;
    int nb_dist = boost::python::len(dist_list);
    stat_tool::wrap_util::auto_ptr_array<const Distribution *> 
      dists(new const Distribution*[nb_dist]);

    for (int i = 0; i < nb_dist; i++)
      dists[i] = extract<const Distribution*>(dist_list[i]);

    const Distribution** d = dists.get();

    MultiPlotSet* ret = p.get_plotable_dists(error, nb_dist, d);
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }


};




// Boost.Python Wrapper export function
void class_distribution()
{

  
  // Distribution base class
  class_< Distribution>("_Distribution")
    .def(self_ns::str(self)) // __str__
    .def( self == self )
    .def( self != self )
    ;


  // Parametric base class
  class_< Parametric, bases< Distribution > >
    ("_Parametric", init< optional< int, int, int, int, double, double > >())
    .def(init<int, int>())
    .def(init<Distribution&>())
    .def(init<Parametric&>())
    .def_readonly("ident", &Parametric::ident)
    .def_readonly("inf_bound", &Parametric::inf_bound)
    .def_readonly("sup_bound", &Parametric::sup_bound)
    .def_readonly("parameter", &Parametric::parameter)
    .def_readonly("probability", &Parametric::probability)
    .def(self_ns::str(self))

    .def("simulate", &Parametric::simulation, 
	 "Simulation one value")
    ;


  enum_<stat_tool::wrap_util::UniqueInt<5, 0> >("DistributionIdentifier")
    .value("NON_PARAMETRIC", NONPARAMETRIC)
    .value("BINOMIAL",BINOMIAL)
    .value("POISSON",POISSON)
    .value("NEGATIVE_BINOMIAL",NEGATIVE_BINOMIAL)
    .value("UNIFORM",UNIFORM)
    .export_values()
    ;
  

  // _Parametric Model
  class_< Parametric_model, bases< Parametric, STAT_interface > >
    ("_ParametricModel", "Parametric model", init <const Histogram& >())

    .def(init< int, int, int, double, double, optional< double > >())
    .def("__init__", make_constructor(ParametricModelWrap::parametric_model_from_file))
    .def(self_ns::str(self)) // __str__ 

    // Output
    .def("get_plotable", ParametricModelWrap::get_plotable,
	 return_value_policy< manage_new_object >(),
	 "Return a plotable for a list of distribution")

    .def("get_plotable", &STAT_interface::get_plotable,
	 return_value_policy< manage_new_object >(),
	 "Return a plotable (no parameters)")

//     .def("plot_write", ParametricModelWrap::plot_write,
// 	 python::args("prefix", "title", "dists"),
// 	 "Write GNUPLOT files (with prefix) for a list of distribution")

//     .def("plot_write", &StatInterfaceWrap::plot_write,
// 	 python::args("prefix", "title"),
// 	  "Write GNUPLOT files (with prefix)")


    .def("survival_ascii_write", ParametricModelWrap::survival_ascii_write,
	 "Return a string containing the object description (survival viewpoint)")

    .def("survival_plot_write", ParametricModelWrap::survival_plot_write,
	 python::args("prefix", "title"),
	 "Write GNUPLOT files (survival viewpoint)")

    .def("survival_get_plotable", ParametricModelWrap::survival_get_plotable,
	  return_value_policy< manage_new_object >(),
	 "Return a plotable object")
    
    .def("survival_spreadsheet_write", &ParametricModelWrap::survival_spreadsheet_write,
	 python::arg("filename"),
	 "Write object to filename (spreadsheet format)")
    
    // Extract
    .def("extract_data", ParametricModelWrap::extract_data, 
	 return_value_policy< manage_new_object >(),
	 "Return the 'data' part of the model")

    .def("simulate", ParametricModelWrap::simulation, 
	 return_value_policy< manage_new_object >(),
	 python::arg("nb_value"),
	 "Simulate values")

    .def("simulate", &Parametric::simulation, 
	 "Simulation one value")

    ;

}
