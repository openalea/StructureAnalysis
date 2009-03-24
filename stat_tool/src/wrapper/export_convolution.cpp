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

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;




////////////////////// Export Convolution ////////////////////////////////////////

class ConvolutionWrap
{

public:

  static boost::shared_ptr<Convolution> convolution_from_file(char* filename)
  {
    Format_error error;
    Convolution *conv = NULL;
    conv = convolution_ascii_read(error, filename);

    if(!conv)
    {
	  stat_tool::wrap_util::throw_error(error);
    }

    return boost::shared_ptr<Convolution>(conv);
  }

  static boost::shared_ptr<Convolution> convolution_from_dists(boost::python::list& dists)
  {
    Format_error error;
    Convolution *conv = NULL;
    int nb_dist = 0;

    nb_dist = boost::python::len(dists);

    if(nb_dist == 0)
    {
      stat_tool::wrap_util::throw_error("Input list cannot be empty");
    }

    stat_tool::wrap_util::auto_ptr_array<const Parametric *>
      dist(new const Parametric*[nb_dist]);

    for(int i=0; i<nb_dist; i++)
	dist[i] = boost::python::extract< Parametric *>(dists[i]);

    conv = convolution_building(error, nb_dist, dist.get());

    if(!conv)
    	stat_tool::wrap_util::throw_error(error);


    return boost::shared_ptr<Convolution>(conv);
  }



  static Convolution_data* simulation(const Convolution& convol, int nb_element)
  {
    Format_error error;
    Convolution_data* ret = NULL;

    ret = convol.simulation(error, nb_element);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Parametric_model* extract(const Convolution& convol, int index)
  {
    Format_error error;
    Parametric_model* ret = NULL;

    ret = convol.extract(error, index);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Parametric_model* extract_convolution(const Convolution& convol)
  {
    Parametric_model* ret;
    Convolution_data* convol_histo = NULL;

    convol_histo = convol.get_convolution_data();
    ret = new Parametric_model(convol,
			       (convol_histo ? convol_histo->get_convolution() : NULL));
    return ret;
  }


  static Convolution_data* extract_data(const Convolution& convol)
  {
    Format_error error;
    Convolution_data* ret = NULL;

    ret = convol.extract_data(error);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static void file_ascii_write(const Convolution& m, const char* path, bool exhaustive)
  {
     bool result = true;
     Format_error error;

     result = m.ascii_write(error, path, exhaustive);
     if (!result)
        stat_tool::wrap_util::throw_error(error);

   }
};



// Boost declaration

void class_convolution()
{

  class_< Convolution, bases< Distribution, STAT_interface > >
    ("_Convolution", "Convolution Distribution")
    .def("__init__", make_constructor(ConvolutionWrap::convolution_from_dists))
    .def("__init__", make_constructor(ConvolutionWrap::convolution_from_file))

    .def(self_ns::str(self)) // __str__

    .def("simulate",   ConvolutionWrap::simulation,
	 return_value_policy< manage_new_object >(),
	 python::arg("nb_element"),
	 "Simulate elements"
	 )

    .def("nb_distribution", &Convolution::get_nb_distribution,
	 "Return the number of components")

    .def("extract_elementary", ConvolutionWrap::extract,
	 return_value_policy< manage_new_object >(),
	 python::arg("index"),
	 "Extract a particular element. First index is 1")

    .def("extract_convolution", ConvolutionWrap::extract_convolution,
	 return_value_policy< manage_new_object >(),
	 "Return a _ParametricModel object")

    .def("extract_data", ConvolutionWrap::extract_data,
	 return_value_policy< manage_new_object >(),
	 "Return the associated _ConvolutionData")

	.def("file_ascii_write", ConvolutionWrap::file_ascii_write,
     "Save Convolution into a file")

    ;
}



////////////////////////// Class Convolution_data //////////////////////////////////

class ConvolutionDataWrap
{

public:

  static Distribution_data* extract(const Convolution_data& convol, int index)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = convol.extract(error, index+1);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Distribution_data* extract_convolution(const Convolution_data& convol_histo)
  {
    Distribution_data* ret;

    ret = new Distribution_data(convol_histo, convol_histo.get_convolution());

    return ret;
  }

};




void class_convolution_data()
{
  class_< Convolution_data, bases< Histogram, STAT_interface > >
    ("_ConvolutionData", "Convolution Data")

    .def(self_ns::str(self))

    .def("nb_histogram", &Convolution_data::get_nb_histogram)
    .def("extract_elementary", ConvolutionDataWrap::extract,
	 return_value_policy< manage_new_object >(),
	 python::arg("index"),
	 "Extract a particular element. First index is 1"
	 )
    .def("extract_convolution", ConvolutionDataWrap::extract_convolution,
	 return_value_policy< manage_new_object >(),
	 "Return a _DistributionData"
	 )

    ;
}




