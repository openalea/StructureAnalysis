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
 *        $Id$
 *
 *----------------------------------------------------------------------------*/



#include "wrapper_util.h"

#include "stat_tool/compound.h"
#include "stat_tool/convolution.h"
#include "stat_tool/discrete_mixture.h"
// #include "stat_tool/multivariate_mixture.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;



////////////////////// Export DiscreteMixture ////////////////////////////////////////
#define WRAP DiscreteMixtureWrap
class WRAP
{

public:

  static boost::shared_ptr<DiscreteMixture>
  mixture_from_file(char *filename)
  {
    StatError error;
    DiscreteMixture *mix = NULL;
    mix = DiscreteMixture::ascii_read(error, filename);
    if (!mix)
      stat_tool::wrap_util::throw_error(error);
    return boost::shared_ptr<DiscreteMixture>(mix);
  }

  static boost::shared_ptr<DiscreteMixture>
  mixture_from_dists(boost::python::list& weights, boost::python::list& dists)
  {
    StatError error;
    DiscreteMixture *mix = NULL;
    int nb_component = 0;

    nb_component = boost::python::len(weights);

    // Test list length
    if (nb_component != boost::python::len(dists))
      {
        stat_tool::wrap_util::throw_error(
            "Input lists must have the same length");
      }
    // Test list length
    if (nb_component == 0)
      {
        stat_tool::wrap_util::throw_error("Input lists cannot be empty");
      }

    stat_tool::wrap_util::auto_ptr_array<double> weight(
        new double[nb_component]);

    stat_tool::wrap_util::auto_ptr_array<const DiscreteParametric*> component(
        new const DiscreteParametric*[nb_component]);

    int i = 0;
    for (i = 0; i < nb_component; i++)
      {
        weight[i] = boost::python::extract<double>(weights[i]);

        boost::python::extract<DiscreteParametric*> get_param(dists[i]);
        if (get_param.check())
          {
            component[i] = new DiscreteParametric(*get_param());
          }
        else
          {
            component[i] = new DiscreteParametric(
                *boost::python::extract<Distribution*>(dists[i])());
          }

      }

    mix = DiscreteMixture::building(error, nb_component, weight.get(), component.get());

    if (!mix)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<DiscreteMixture>(mix);
  }

  static boost::shared_ptr<DiscreteMixture>
  mixture_from_unknown_component(boost::python::list& dists)
  {
    StatError error;
    DiscreteMixture *mix = NULL;
    int nb_component = 0;

    nb_component = boost::python::len(dists);

    // Test list length
    if (nb_component == 0)
      stat_tool::wrap_util::throw_error("Input list cannot be empty");

    stat_tool::wrap_util::auto_ptr_array<const DiscreteParametric *> component(
        new const DiscreteParametric*[nb_component]);

    for (int i = 0; i < nb_component; i++)
      {
        component[i] = boost::python::extract<DiscreteParametric *>(dists[i]);
      }

    mix = new DiscreteMixture(nb_component, component.get());
    if (!mix)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<DiscreteMixture>(mix);
  }

  static DiscreteParametricModel*
  extract_weight(const DiscreteMixture &mixt)
  {
    StatError error;
    DiscreteParametricModel *ret;
    DiscreteMixtureData *mixt_data = NULL;

    mixt_data = mixt.get_mixture_data();
    ret = new DiscreteParametricModel(*(mixt.get_weight()),
        (mixt_data ? mixt_data->get_weight() : NULL));
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static DiscreteParametricModel*
  extract_mixture(DiscreteMixture &mixture_input)
  {
    StatError error;
    DiscreteParametricModel *ret;
    DiscreteMixtureData *mixture_data = NULL;

    mixture_data = mixture_input.get_mixture_data();

    ret = new DiscreteParametricModel(*((Distribution*) (&mixture_input)),
        (FrequencyDistribution*) mixture_data);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static DiscreteParametricModel*
  extract_component(const DiscreteMixture &input, int var1)
  {
    StatError error;
    DiscreteParametricModel *ret = NULL;

    ret = input.extract(error, var1);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // component case
  static DiscreteMixtureData*
  extract_data(const DiscreteMixture &input)
  {
    StatError error;
    DiscreteMixtureData* ret = NULL;

    ret = input.extract_data(error);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  WRAP_METHOD1(DiscreteMixture, simulation, DiscreteMixtureData, int);
  WRAP_METHOD_FILE_ASCII_WRITE( DiscreteMixture);
  WRAP_METHOD_SPREADSHEET_WRITE( DiscreteMixture);

  static MultiPlotSet*
  get_plotable(const DiscreteMixture &mixt)
  {
    StatError error;
    MultiPlotSet *ret = mixt.get_plotable();
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static MultiPlotSet*
  survival_get_plotable(const DiscreteMixture &p)
  {
    StatError error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if (!ret)
      ERROR;
    return ret;
  }

};



// Boost declaration
void class_mixture()
{

  class_< DiscreteMixture, bases< Distribution, StatInterface > >
  ("_DiscreteMixture", "DiscreteMixture Distribution")
  // constructors

  .def("__init__", make_constructor(DiscreteMixtureWrap::mixture_from_file),
      "Build from a filename" )
  .def("__init__", make_constructor(DiscreteMixtureWrap::mixture_from_dists),
      "Build from a list of weights and a list of distributions")
  .def("__init__", make_constructor(DiscreteMixtureWrap::mixture_from_unknown_component),
      "Build from unknown components") // internal use

  // Python Operators
  .def("__len__", &DiscreteMixture::get_nb_component,
      "Return the number of components") // __len__
  .def(self_ns::str(self)) // __str__

  // properties
  .add_property("nb_component", &DiscreteMixture::get_nb_component,
      "Return the number of components")

  // python modules
  DEF_RETURN_VALUE("simulate", WRAP::simulation,
      args("nb_element"), "See Simulate")

  // Used in Python modules such as:
  // ----------------------- Extract
  DEF_RETURN_VALUE("extract_component", WRAP::extract_component,
      args("index"), "Get a particular component. First index is 1")
  DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight,
      "Return the weight distribution")
  DEF_RETURN_VALUE_NO_ARGS("extract_mixture", WRAP::extract_mixture,
      "Return the DiscreteMixture distribution")
  DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,
      "Return the associated _DiscreteMixtureData object")


  //others
  //  DEF_RETURN_VALUE_NO_ARGS("get_plotable", &StatInterface::get_plotable,"Return a plotable (no parameters)");
  .def("file_ascii_write", WRAP::file_ascii_write,
      "Save Mixture into a file")
  .def("spreadsheet_write", WRAP::spreadsheet_write,
      "save data in spreadsheet format")
  DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "return plotable")
  DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable,
      "Return a survival plotable")

  ;

  /*

  DiscreteMixture();
  DiscreteMixture(const DiscreteMixture &mixt , bool *component_flag , int inb_value);
  DiscreteMixture(int inb_component , const DiscreteParametric **pcomponent);
  DiscreteMixture(const DiscreteMixture &mixt , bool data_flag = true)    :Distribution(mixt)
    { copy(mixt , data_flag); }


  void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD ,
                     bool component_flag = true);
  double likelihood_computation(const DiscreteMixtureData &mixt_histo) const;

  // acces membres de la classe

  DiscreteParametric* get_component(int index) const { return component[index]; }
*/
}
#undef WRAP


////////////////////////// Class DiscreteMixtureData //////////////////////////////////
#define WRAP DiscreteMixtureDataWrap
class DiscreteMixtureDataWrap
{

public:

  WRAP_METHOD1(DiscreteMixtureData, extract, DiscreteDistributionData, int);
  WRAP_METHOD_FILE_ASCII_WRITE( DiscreteMixtureData);

  static DiscreteDistributionData*
  extract_weight(const DiscreteMixtureData &mixt_histo)
  {
    DiscreteDistributionData *ret;
    ret = new DiscreteDistributionData(*(mixt_histo.get_weight()),
        mixt_histo.get_mixture()->get_weight());
    return ret;
  }

  static DiscreteDistributionData*
  extract_mixture(const DiscreteMixtureData &mixt_histo)
  {
    DiscreteDistributionData *ret;
    ret = new DiscreteDistributionData(mixt_histo, mixt_histo.get_mixture());
    return ret;
  }

};

void class_mixture_data()
{
  class_< DiscreteMixtureData, bases< FrequencyDistribution, StatInterface > >
  ("_DiscreteMixtureData",  "DiscreteMixture Data")

  // Python Operators
  .def(self_ns::str(self)) //str

  // properties
 .add_property("nb_component", &DiscreteMixtureData::get_nb_component,
      "Return the number of components.")

  // getters
  DEF_RETURN_VALUE("get_component", &DiscreteMixtureData::get_component,
      args("index"), "Return the number of components.")

  // Used in Python modules such as:
  // -----------------------  Extract
  DEF_RETURN_VALUE("extract_component", WRAP::extract,
      args("index"),"Get a particular component. First index is 1")
  DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight,
      "Return a _DistributionData")
  DEF_RETURN_VALUE_NO_ARGS("extract_mixture", WRAP::extract_mixture,
      "Return a _DistributionData")


  //others

  .def("file_ascii_write", WRAP::file_ascii_write, "Save Compound into a file")

  /*
    DiscreteMixtureData(const FrequencyDistribution &histo , int inb_component);
    DiscreteMixtureData(const DiscreteMixture &mixt);
    DiscreteMixtureData(const DiscreteMixtureData &mixt_histo , bool model_flag = true) :FrequencyDistribution(mixt_histo) { copy(mixt_histo , model_flag); }

    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix , const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    FrequencyDistribution* get_component(int index) const { return component[index]; }
    */
    ;
}

#undef WRAP


