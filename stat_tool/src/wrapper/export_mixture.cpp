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

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/convolution.h"
#include "stat_tool/mixture.h"
#include "stat_tool/compound.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;



////////////////////// Export Mixture ////////////////////////////////////////

class MixtureWrap
{

public:

  static boost::shared_ptr<Mixture> mixture_from_file(char* filename)
  {
    Format_error error;
    Mixture *mix = NULL;
    mix = mixture_ascii_read(error, filename);

    if(!mix)
      {
	stat_tool::wrap_util::throw_error(error);
      }

    return boost::shared_ptr<Mixture>(mix);
  }

  static boost::shared_ptr<Mixture> mixture_from_components(boost::python::list& weights, 
							    boost::python::list& dists)
  {
    Format_error error;
    Mixture *mix = NULL;
    int nb_component = 0; 
    
    nb_component = boost::python::len(weights);
    
    // Test list length
    if(nb_component != boost::python::len(dists))
      {
	stat_tool::wrap_util::throw_error("Input lists must have the same length");
      }
    // Test list length
    if(nb_component == 0)
      {
	stat_tool::wrap_util::throw_error("Input lists cannot be empty");
      }

    
    double* weight = new double[nb_component];
    const Parametric **component = new const Parametric*[nb_component];

    try
      {
	for(int i=0; i<nb_component; i++)
	  {
	    weight[i] = boost::python::extract< double >(weights[i]);
	    component[i] = boost::python::extract< Parametric *>(dists[i]);
	    
	  }
      }
    catch(...)
      {
	delete[] weight;
	delete[] component;
	throw;
      }
    
    
    mix = mixture_building(error, nb_component, weight, component);
    
    delete[] weight;
    delete[] component;
	
    if(!mix)
      {
	stat_tool::wrap_util::throw_error(error);
      }
    
    return boost::shared_ptr<Mixture>(mix);
  }


  static boost::shared_ptr<Mixture> mixture_from_unknown_component(boost::python::list& dists)
  {
    Format_error error;
    Mixture *mix = NULL;
    int nb_component = 0; 
    
    nb_component = boost::python::len(dists);
    
    // Test list length
    if(nb_component == 0)
      {
	stat_tool::wrap_util::throw_error("Input list cannot be empty");
      }

    
    const Parametric **component = new const Parametric*[nb_component];

    try
      {
	for(int i=0; i<nb_component; i++)
	  component[i] = boost::python::extract< Parametric *>(dists[i]);
      }
    catch(...)
      {
	delete[] component;
	throw;
      }
    
    
    mix = new Mixture(nb_component, component);
    
    delete[] component;
	
    return boost::shared_ptr<Mixture>(mix);
  }

  
  static Mixture_data* simulation(const Mixture& mixt, int nb_element)
  {
    Format_error error;
    Mixture_data* ret = NULL;

    ret = mixt.simulation(error, nb_element);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Parametric_model* extract(const Mixture& mixt, int index)
  {
    Format_error error;
    Parametric_model* ret = NULL;

    ret = mixt.extract(error, index);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Parametric_model* extract_weight(const Mixture& mixt)
  {
    Parametric_model* ret;
    Mixture_data* mixt_histo = NULL;

    mixt_histo = mixt.get_mixture_data();
    ret = new Parametric_model(*(mixt.get_weight()),
			       (mixt_histo ? mixt_histo->get_weight() : NULL));
    return ret;
  }


  static Parametric_model* extract_mixture(const Mixture& mixt)
  {
    Parametric_model* ret;
    Mixture_data* mixt_histo = NULL;

    mixt_histo = mixt.get_mixture_data();
    ret = new Parametric_model(mixt,
			       (mixt_histo ? mixt_histo->get_mixture() : NULL));
    return ret;
  }


  static Mixture_data* extract_data(const Mixture& mixt)
  {
    Format_error error;
    Mixture_data* ret = NULL;

    ret = mixt.extract_data(error);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


};



// Boost declaration

void class_mixture()
{

  class_< Mixture, bases< Distribution, STAT_interface > >
    ("_Mixture", "Mixture Distribution")
    .def("__init__", make_constructor(MixtureWrap::mixture_from_file),
	 "Build from a filename"
	 )
    
    .def("__init__", make_constructor(MixtureWrap::mixture_from_components),
	 "Build from a list of weights and a list of distributions")
	 

    .def("__init__", make_constructor(MixtureWrap::mixture_from_unknown_component),
	 "Build from unknown components") // internal use

    .def(self_ns::str(self)) // __str__ 
    
    .def("simulate", MixtureWrap::simulation, 
	 return_value_policy< manage_new_object >(),
	 python::arg("nb_element"),
	 "Simulate nb_element elements")

    .def("nb_component", &Mixture::get_nb_component,
	 "Return the number of components")
    
    .def("extract_component", MixtureWrap::extract, 
	 return_value_policy< manage_new_object >(),
	 python::arg("index"),
	 "Get a particular component. First index is 1")

    .def("extract_weight", MixtureWrap::extract_weight, 
	 return_value_policy< manage_new_object >(),
	 "Return the weight distribution")

    .def("extract_mixture", MixtureWrap::extract_mixture, 
	 return_value_policy< manage_new_object >(),
	 "Return the Mixture distribution")

    .def("extract_data", MixtureWrap::extract_data, 
	 return_value_policy< manage_new_object >(),
	 "Return the associated _MixtureData object"
	 )
    ;
}



////////////////////////// Class Mixture_data //////////////////////////////////

class MixtureDataWrap
{

public:

  static Distribution_data* extract(const Mixture_data& d, int index)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = d.extract(error, index);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data* extract_weight(const Mixture_data& mixt_histo)
  {
    Distribution_data* ret;

    ret = new Distribution_data(*(mixt_histo.get_weight()) ,
				mixt_histo.get_mixture()->get_weight());

    return ret;
  }


  static Distribution_data* extract_mixture(const Mixture_data& mixt_histo)
  {
    Distribution_data* ret;

    ret = new Distribution_data(mixt_histo , mixt_histo.get_mixture());

    return ret;
  }


};




void class_mixture_data()
{
  class_< Mixture_data, bases< Histogram, STAT_interface > >
    ("_MixtureData",  "Mixture Data")

    .def(self_ns::str(self))

    .def("nb_component", &Mixture_data::get_nb_component,
	 "Return the number of components."
	 )

    .def("extract_component", MixtureDataWrap::extract, 
	 return_value_policy< manage_new_object >(),
	 python::arg("index"),
	 "Get a particular component. First index is 1"
	 )
    
    .def("extract_weight", MixtureDataWrap::extract_weight, 
	 return_value_policy< manage_new_object >(),
	 "Return a _DistributionData"
	 )
    
    .def("extract_mixture", MixtureDataWrap::extract_mixture, 
	 return_value_policy< manage_new_object >(),
	 "Return a _DistributionData"
	 )
    
    ;
}




