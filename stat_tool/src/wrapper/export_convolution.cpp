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
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id: export_convolution.cpp 18031 2015-04-23 07:13:15Z guedon $
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
using namespace stat_tool;



////////////////////// Export Convolution ////////////////////////////////////////

class ConvolutionWrap
{

public:

  static boost::shared_ptr<Convolution>
  convolution_from_file(char* filename)
  {
    StatError error;
    Convolution *conv = NULL;
    conv = convolution_ascii_read(error, filename);

    if (!conv)
      {
        stat_tool::wrap_util::throw_error(error);
      }

    return boost::shared_ptr<Convolution>(conv);
  }

  static boost::shared_ptr<Convolution>
  convolution_from_dists(boost::python::list& dists)
  {
    StatError error;
    Convolution *conv = NULL;
    int nb_dist = 0;

    nb_dist = boost::python::len(dists);
    if (nb_dist == 0)
      {
        stat_tool::wrap_util::throw_error("Input list cannot be empty");
      }

    stat_tool::wrap_util::auto_ptr_array<const DiscreteParametric *> dist(
        new const DiscreteParametric*[nb_dist]);

    int i = 0;
    for (i = 0; i < nb_dist; i++)
      {

        boost::python::extract<DiscreteParametric*> get_param(dists[i]);
        if (get_param.check())
          {
            dist[i] = new DiscreteParametric(*get_param());
          }
        else
          {
            dist[i] = new DiscreteParametric(*boost::python::extract<Distribution*>(
                dists[i])());
          }
      }

    conv = convolution_building(error, nb_dist, dist.get());

    if (!conv)
      stat_tool::wrap_util::throw_error(error);

    for (i = 0; i < nb_dist; i++)
      {
        delete dist[i];
      }

    return boost::shared_ptr<Convolution>(conv);
  }

  WRAP_METHOD1(Convolution, simulation, ConvolutionData, int); // simulate
  WRAP_METHOD1(Convolution, extract, DiscreteParametricModel, int); //extract_elementary
  WRAP_METHOD0(Convolution, extract_data, ConvolutionData); //extract_data
  WRAP_METHOD_FILE_ASCII_WRITE( Convolution);

  static DiscreteParametricModel*
  extract_elementary(const Convolution& input, int index)
  {
    StatError error;
    DiscreteParametricModel *ret = NULL;
    ret = input.extract(error, index);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static DiscreteParametricModel*
  extract_convolution(const Convolution& convolution_input)
  {
    DiscreteParametricModel *ret;
    ConvolutionData *convolution_data = NULL;
    convolution_data = convolution_input.get_convolution_data();
    ret = new DiscreteParametricModel(*((Distribution*) (&convolution_input)),
        (FrequencyDistribution*) convolution_data);
    return ret;
  }

  static MultiPlotSet*
  survival_get_plotable(const Convolution &p)
  {
    StatError error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if (!ret)
      ERROR;
    return ret;
  }

  static MultiPlotSet*
  get_plotable(const Convolution &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }

};



// Boost declaration

#define WRAP ConvolutionWrap
void class_convolution()
{

  class_< Convolution, bases< Distribution, StatInterface > >
    ("_Convolution", "Convolution Distribution")
    .def("__init__", make_constructor(WRAP::convolution_from_dists))
    .def("__init__", make_constructor(WRAP::convolution_from_file))
    .def("__len__", &Convolution::get_nb_distribution,
        "Return the size of the Class instance")
    .def(self_ns::str(self))
    .def("nb_distribution", &Convolution::get_nb_distribution,
        "Return the number of components")

    DEF_RETURN_VALUE("simulate", WRAP::simulation,
        args("nb_element"), "Simulate elements")
    // check extract and extract_data
    DEF_RETURN_VALUE("extract", WRAP::extract,
        args("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE("extract_elementary", WRAP::extract_elementary, args("index"),
        "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE_NO_ARGS("extract_convolution", WRAP::extract_convolution,
        "Return a _ParametricModel object")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,
        "Return the associated _ConvolutionData")
    DEF_RETURN_VALUE_NO_ARGS("file_ascii_write", WRAP::file_ascii_write,
        "Save Convolution into a file")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable,
        "Return a survival plotable")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, ""
        "Return a plotable")
    ;
/*

   Convolution(const Convolution &convol , bool data_flag = true):Distribution(convol) { copy(convol , data_flag); }
   void computation(int min_nb_value = 1 , double cumul_threshold = CONVOLUTION_THRESHOLD ,  bool *dist_flag = 0);
   DiscreteParametric* get_distribution(int index) const { return distribution[index]; }
   */



}
#undef WRAP


#define WRAP ConvolutionDataWrap

////////////////////////// Class ConvolutionData //////////////////////////////////

class ConvolutionDataWrap
{

public:

  static DiscreteDistributionData* extract(const ConvolutionData &convol_histo, int index)
  {
    StatError error;
    DiscreteDistributionData* ret = NULL;

    ret = convol_histo.extract(error, index);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static DiscreteDistributionData* extract_convolution(const ConvolutionData &convol_histo)
  {
    DiscreteDistributionData* ret;
    ret = new DiscreteDistributionData(convol_histo, convol_histo.get_convolution());
    return ret;
  }

};




void class_convolution_data()
{
  class_< ConvolutionData, bases< FrequencyDistribution, StatInterface > >
    ("_ConvolutionData", "Convolution Data")
    .def(self_ns::str(self))
    .def("nb_distribution", &ConvolutionData::get_nb_distribution)
    .def("get_frequency_distribution", &ConvolutionData::get_frequency_distribution,
        return_value_policy< manage_new_object >(), args("index"),"todo")
    DEF_RETURN_VALUE("extract", ConvolutionDataWrap::extract,
        args("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE("extract_elementary", ConvolutionDataWrap::extract,
        args("index"), "Extract a particular element. First index is 1")
    DEF_RETURN_VALUE_NO_ARGS("extract_convolution", ConvolutionDataWrap::extract_convolution,
        "Return a _DistributionData")
    ;

  /*
    ConvolutionData(const FrequencyDistribution &histo , int nb_histo);
    ConvolutionData(const Convolution &convol);
    ConvolutionData(const ConvolutionData &convol_histo , bool model_flag = true):FrequencyDistribution(convol_histo) { copy(convol_histo , model_flag); }
    FrequencyDistribution* get_frequency_distribution(int index) const { return histogram[index]; }
    */



}

#undef WRAP


