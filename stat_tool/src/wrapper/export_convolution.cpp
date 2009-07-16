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


#include "boost_python_aliases.h"

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
    if(nb_dist == 0){
      stat_tool::wrap_util::throw_error("Input list cannot be empty");
    }

    stat_tool::wrap_util::auto_ptr_array<const Parametric *>
      dist(new const Parametric*[nb_dist]);

    for(int i=0; i<nb_dist; i++){
    	dist[i] = boost::python::extract< Parametric *>(dists[i]);
    }

    conv = convolution_building(error, nb_dist, dist.get());

    if(!conv)
    	stat_tool::wrap_util::throw_error(error);


    return boost::shared_ptr<Convolution>(conv);
  }



 WRAP_METHOD1(Convolution, simulation, Convolution_data, int);  // simulate
 WRAP_METHOD1(Convolution, extract, Parametric_model, int); 	//extract_elementary
 WRAP_METHOD0(Convolution, extract_data, Convolution_data);		//extract_data
 WRAP_METHOD_FILE_ASCII_WRITE(Convolution);



  static Parametric_model* extract_elementary(const Convolution& input, int index) 
  {
    Format_error error;
    Parametric_model* ret = NULL;
    
    ret = input.extract(error, index); 
    if(!ret) 
        stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }
  
  static Parametric_model* extract_convolution(const Convolution& convolution_input)
  {
    Parametric_model* ret;
    Convolution_data* convolution_data = NULL;
    
    convolution_data = convolution_input.get_convolution_data();
    
    //ret = new Parametric_model(convolution_input,
	//		       (convolution_data ? convolution_data->get_convolution() : NULL));
			       
    ret = new Parametric_model(*((Distribution*)(&convolution_input)), 
            (Histogram*)convolution_data);			       
    		       
    return ret;
  }

  static MultiPlotSet* survival_get_plotable(const Convolution& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if (!ret) ERROR;
    return ret;
  }
 
  static MultiPlotSet* get_plotable(const Convolution& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret) ERROR;
    return ret;
  }

};



// Boost declaration

#define WRAP ConvolutionWrap
void class_convolution()
{

  class_< Convolution, bases< Distribution, STAT_interface > >
    ("_Convolution", "Convolution Distribution")
    .def("__init__", make_constructor(WRAP::convolution_from_dists))
    .def("__init__", make_constructor(WRAP::convolution_from_file))
    .def("__len__", &Convolution::get_nb_distribution, "Return the size of the Class instance")
    .def(self_ns::str(self))
    .def("nb_distribution", &Convolution::get_nb_distribution, "Return the number of components")

	DEF_RETURN_VALUE("simulate", WRAP::simulation, ARGS("nb_element"), "Simulate elements")
	// check extract and extract_data 
    DEF_RETURN_VALUE("extract", WRAP::extract, ARGS("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE("extract_elementary", WRAP::extract_elementary, ARGS("index"), "Extract a particular element. First index is 1")
	DEF_RETURN_VALUE_NO_ARGS("extract_convolution", WRAP::extract_convolution, "Return a _ParametricModel object")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data, "Return the associated _ConvolutionData")
    DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write, "Save Convolution into a file")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable, "Return a survival plotable")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")
    ;
/*

   Convolution(const Convolution &convol , bool data_flag = true):Distribution(convol) { copy(convol , data_flag); }

   void computation(int min_nb_value = 1 , double cumul_threshold = CONVOLUTION_THRESHOLD ,  bool *dist_flag = 0);
   Convolution_data* simulation(Format_error &error , int nb_element) const;

   Parametric* get_distribution(int index) const { return distribution[index]; }
   */



}
#undef WRAP


#define WRAP ConvolutionDataWrap

////////////////////////// Class Convolution_data //////////////////////////////////

class ConvolutionDataWrap
{

public:

  static Distribution_data* extract(const Convolution_data& convol, int index)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = convol.extract(error, index);
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
    .def("get_histogram", &Convolution_data::get_histogram,	return_value_policy< manage_new_object >(), ARGS("index"),"todo")
    DEF_RETURN_VALUE("extract", ConvolutionDataWrap::extract, ARGS("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE("extract_elementary", ConvolutionDataWrap::extract, ARGS("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE_NO_ARGS("extract_convolution", ConvolutionDataWrap::extract_convolution,"Return a _DistributionData")
    ;

  /*
    Convolution_data(const Histogram &histo , int nb_histo);
    Convolution_data(const Convolution &convol);
    Convolution_data(const Convolution_data &convol_histo , bool model_flag = true)
    :Histogram(convol_histo) { copy(convol_histo , model_flag); }


    Histogram* get_histogram(int index) const { return histogram[index]; }
    */



}

#undef WRAP


