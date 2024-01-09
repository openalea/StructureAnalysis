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
 *        $Id: export_mixture.cpp 18037 2015-04-23 07:16:12Z guedon $
 *
 *----------------------------------------------------------------------------*/



#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/compound.h"
#include "stat_tool/convolution.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/multivariate_mixture.h"

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
    mix = discrete_mixture_ascii_read(error, filename);
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

    mix = discrete_mixture_building(error, nb_component, weight.get(), component.get());

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



////////////////////////// Class MultivariateMixture //////////////////////////////////

#define WRAP MultivariateMixtureWrap
class MultivariateMixtureWrap
{

public:

  static boost::shared_ptr<MultivariateMixture>
  multivariate_mixture_from_file(char* filename)
  {
    StatError error;
    MultivariateMixture *mix = NULL;
    mix = multivariate_mixture_ascii_read(error, filename);

    if (mix == NULL)
      {
        stat_tool::wrap_util::throw_error(error);
      }

    return boost::shared_ptr<MultivariateMixture>(mix);
  }

  static boost::shared_ptr<MultivariateMixture>
  multivariate_mixture_from_mixture(const MultivariateMixture& mixt)
  {
    MultivariateMixture *mix_cp = NULL;

    mix_cp = new MultivariateMixture(mixt, true);

    return boost::shared_ptr<MultivariateMixture>(mix_cp);
  }

  static boost::shared_ptr<MultivariateMixture>
  multivariate_mixture_from_components(boost::python::list& weights,
      boost::python::list& dists)
  {
    StatError error;
    MultivariateMixture *mix = NULL;
    int nb_component = 0;
    int nb_variable = 0;
    // DiscreteParametric *pcomp = NULL;
    // boost::python::list comp_list;


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

    nb_variable = boost::python::len(dists[0]);

    if (nb_variable == 0)
      {
        stat_tool::wrap_util::throw_error("Number of variable must be positive");
      }

    stat_tool::wrap_util::auto_ptr_array<double> weight(
        new double[nb_component]);

    stat_tool::wrap_util::auto_ptr_array<DiscreteParametricProcess *> pcomponent(
        new DiscreteParametricProcess*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<CategoricalProcess *> npcomponent(
        new CategoricalProcess*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<bool> is_parametric(
        new bool[nb_variable]);

    int i, var;

    for (i = 0; i < nb_component; i++)
      weight[i] = boost::python::extract<double>(weights[i]);

    for (var = 0; var < nb_variable; var++)
      {

        DiscreteParametric **pprocess = new DiscreteParametric*[nb_component];
        Distribution **npprocess = new Distribution*[nb_component];

        pprocess[0] = NULL;
        npprocess[0] = NULL;

        for (i = 0; i < nb_component; i++)
          {
            boost::python::extract<DiscreteParametric*> x(dists[i][var]);
            if (x.check())
              {
                pprocess[i] = x();
                npprocess[i] = NULL;
                if (i == 0)
                  is_parametric[var] = true;
                else
                  {
                    if (!is_parametric[var])
                      stat_tool::wrap_util::throw_error(
                          "Distributions must be the same type for a given variable");
                  }
              }
            else
              {
                boost::python::extract<Distribution *> x(dists[i][var]);
                if (x.check())
                  {
                    npprocess[i] = x();
                    pprocess[i] = NULL;
                    if (i == 0)
                      is_parametric[var] = false;
                    else
                      {
                        if (is_parametric[var])
                          stat_tool::wrap_util::throw_error(
                              "Distributions must be the same type for a given variable");
                      }
                  }
                else
                  stat_tool::wrap_util::throw_error(
                      "Bad type for list element: must be DiscreteParametric or Distribution");

              }
          } // end for (i)
        if (pprocess[0] != NULL)
          {
            pcomponent[var] = new DiscreteParametricProcess(nb_component, pprocess);
            npcomponent[var] = NULL;
          }
        else
          {
            pcomponent[var] = NULL;
            npcomponent[var] = new CategoricalProcess(nb_component,
                npprocess);
          }

        /*for(i=0; i<nb_component; i++) {
         if (npprocess[i] != NULL)
         delete npprocess[i];
         if (pprocess[i] != NULL)
         delete pprocess[i];
         }*/

        delete[] pprocess;
        delete[] npprocess;
      } // end for (var)

    mix = multivariate_mixture_building(error, nb_component, nb_variable, weight.get(),
        pcomponent.get(), npcomponent.get());

    for (var = 0; var < nb_variable; var++)
      {
        if (pcomponent[var] != NULL)
          delete pcomponent[var];
        if (npcomponent[var] != NULL)
          delete npcomponent[var];
      }

    if (mix == NULL)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<MultivariateMixture>(mix);
  }

  static DiscreteParametricModel*
  extract_weight(const MultivariateMixture& mixt)
  {
    DiscreteParametricModel* ret;
    MultivariateMixtureData* mixt_data = NULL;

    mixt_data = mixt.get_mixture_data();
    ret = new DiscreteParametricModel(*(mixt.get_weight()),
        (mixt_data ? mixt_data->get_weight() : NULL));
    return ret;
  }

  static DiscreteParametricModel* 
  extract_mixture(const MultivariateMixture& mixt, int ivariable)
  {
    StatError error;
    DiscreteParametricModel* ret;
    MultivariateMixtureData* mixt_data = NULL;
    Distribution *marginal = mixt.extract_distribution(error, ivariable);
    FrequencyDistribution *marginal_hist = NULL;

    if (marginal == NULL)
      stat_tool::wrap_util::throw_error(error);

    mixt_data = mixt.get_mixture_data();

    if (mixt_data != NULL)
      marginal_hist = mixt_data->get_marginal_distribution(ivariable);

    ret = new DiscreteParametricModel(*marginal, marginal_hist);

    // margina_hist is only a reference, included in mixt_data
    /* if (marginal_hist != NULL)
     delete marginal_hist;*/

    if (marginal != NULL)
      {
        delete marginal;
        marginal = NULL;
      }

    return ret;
  }

  static MultivariateMixtureData* 
  cluster_data(const MultivariateMixture& mixt, Vectors& vec, bool state_entropy=false)
  {
    StatError error;
    MultivariateMixtureData* ret = NULL;

    ret = mixt.cluster(error, vec, VITERBI, state_entropy);

    if (ret == NULL)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static bool
  _is_parametric(const MultivariateMixture& mixt, int ivariable)
  {
    ostringstream error_message;

    if ((ivariable < 0) || (ivariable > mixt.get_nb_component()))
      {
        error_message << "Bad variable index: " << ivariable;
        PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
        throw_error_already_set();
      }
    else
      return mixt.is_parametric(ivariable);
  }

  static void
  state_permutation(const MultivariateMixture& mix, boost::python::list perm)
  {
    bool status = true, several_errors = false;
    int llength, i;
    ostringstream error_message;
    object o;
    StatError error;
    int *iperm = NULL;

    llength = boost::python::len(perm);
    if ((llength > 0) && (llength == mix.get_nb_component()))
      {
        iperm = new int[llength];
        for (i = 0; i < llength; i++)
          {
            o = perm[i];
            try
              {
                extract<int> x(o);
                if (x.check())
                  iperm[i] = x();
                else
                  status = false;
              }
            catch (...)
              {
                status = false;
              }
            if (!status)
              {
                if (several_errors)
                  error_message << endl;
                else
                  several_errors = true;
                error_message << "incorrect type for element " << i
                    << " of argument list: expecting an int ";
              }
          }
        if (!status)
          {
            delete[] iperm;
            iperm = NULL;
            PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
            throw_error_already_set();
          }
      }
    else
      {
        status = false;
        error_message << "incorrect permutation" << endl;
        PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
        throw_error_already_set();
      }
    if (status)
      {
        mix.state_permutation(error, iperm);
        if (error.get_nb_error() > 0)
          {
            error_message << error;
            stat_tool::wrap_util::throw_error(error);
          }
        delete[] iperm;
        iperm = NULL;
      }
  }
  

  static MultiPlotSet*
  get_plotable(const MultivariateMixture &mixt)
  {
    StatError error;
    MultiPlotSet *ret = mixt.get_plotable();
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

    WRAP_METHOD0(MultivariateMixture, extract_data, MultivariateMixtureData);
    WRAP_METHOD_FILE_ASCII_WRITE( MultivariateMixture);
    WRAP_METHOD_PLOT_WRITE( MultivariateMixture);
    WRAP_METHOD_SPREADSHEET_WRITE( MultivariateMixture);
    WRAP_METHOD1(MultivariateMixture, simulation, MultivariateMixtureData, int);

};



// Boost declaration



void class_multivariate_mixture()
{
  class_< MultivariateMixture, bases<StatInterface> >
    ("_MultivariateMixture", "Multivariate Mixture Distribution")

    // constructors
    .def("__init__", make_constructor(MultivariateMixtureWrap::multivariate_mixture_from_file),
        "Build from a filename")
    .def("__init__", make_constructor(MultivariateMixtureWrap::multivariate_mixture_from_mixture),
        "Build from a _MultivariateMixture object")
    .def("__init__", make_constructor(MultivariateMixtureWrap::multivariate_mixture_from_components),
        "Build from a list of weights and a list of list of distributions\n" "(components and variables)")

    // Python Operators
    .def(self_ns::str(self)) // __str__
    .def("__len__", &MultivariateMixture::get_nb_component,
        "Return the number of components") // __len__

    // properties
    .add_property("nb_variable", &MultivariateMixture::get_nb_variable,
        "Return the number of variables")
    .add_property("nb_component", &MultivariateMixture::get_nb_component,
        "Return the number of variables")

    // return object and no arguments needed
    DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight,
        "Return the weight distribution")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,
        "Return the associated _MultivariateMixtureData object")


    // return object and arguments needed
    DEF_RETURN_VALUE("simulate", WRAP::simulation,
        args("nb_element"), "simulate(self, nb_element) -> _MultivariateMixtureData. \n\n" "Simulate nb_element elements")
    DEF_RETURN_VALUE("extract_mixture", WRAP::extract_mixture,
        args("variable"), "extract_mixture(self, variable) -> _Distribution. \n\n" "Return the _MultivariateMixture distribution")

    // no object returned, args required
    .def("_is_parametric", WRAP::_is_parametric,
        args("variable"),"_is_parametric(self, variable) -> bool. \n\n" "Return True if the variable is parametric")
    .def("file_ascii_write", WRAP::file_ascii_write,
        args("path", "exhaustive_flag"), "file_ascii_write(self, path, exhaustive_flag) -> None. \n\n""Save _MultivariateMixture into a file")
    .def("plot_write", WRAP::plot_write,
        args("prefix", "title"),"plot_write(self, prefix, title) -> None. \n\n" "Write GNUPLOT files")

    // no object returned, no arguments required
    .def("state_permutation", WRAP::state_permutation, "permutation of the model states")
    .def("spreadsheet_write", WRAP::spreadsheet_write, "save data in spreadsheet format")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "return plotable")

    .def("cluster_data", WRAP::cluster_data, "cluster_data(self, Vectors, entropy) -> _MultivariateMixtureData. \n\n" "Cluster data using the _MultivariateMixture model", return_value_policy< manage_new_object >())


     ;

    /*

    DiscreteParametricModel* extract_parametric_model(StatError &error , int ivariable,    int index) const;
    Distribution* extract_nonparametric_model(StatError &error , int ivariable, int index) const;
    Distribution* extract_distribution(StatError &error , int ivariable) const;
    std::ostream& line_write(std::ostream &os) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,  bool exhaustive = false) const;
    plotable::MultiPlotSet* get_plotable() const;
    double likelihood_computation(const Vectors &mixt_data, bool log_computation=false) const;
    MultivariateMixtureData* cluster(StatError &error,  const Vectors &vec, int algorithm=VITERBI) const;
    bool is_parametric(int ivariable) const;
  MultivariateMixtureData* get_mixture_data() const { return mixture_data; }
  DiscreteParametric* get_weight() const { return weight; }
  DiscreteParametricProcess* get_parametric_process(int variable) const;
  CategoricalProcess* get_categorical_process(int variable) const;
  DiscreteParametric* get_parametric_component(int variable, int index) const;
  Distribution* get_nonparametric_component(int variable, int index) const;
  
 
*/

}

#undef WRAP

////////////////////////// Class MultivariateMixtureData //////////////////////////////////

#define WRAP MultivariateMixtureDataWrap
class WRAP
{

public:

  static boost::shared_ptr<MultivariateMixtureData>
  multivariate_mixture_data_from_mixture_data(const MultivariateMixtureData& mixt)
  {
    MultivariateMixtureData *mix_cp = NULL;

    mix_cp = new MultivariateMixtureData(mixt, true);

    return boost::shared_ptr<MultivariateMixtureData>(mix_cp);
  }

  static DiscreteDistributionData*
  extract_weight(const MultivariateMixtureData& mixt_histo)
  {
    DiscreteDistributionData* ret;

    ret = new DiscreteDistributionData(*(mixt_histo.get_weight()),
        mixt_histo.get_mixture()->get_weight());

    return ret;
  }

  static MultivariateMixture*
  extract_mixture(const MultivariateMixtureData& mixt_histo)
  {
    MultivariateMixture* ret, *cp_ret = NULL;

    cp_ret = mixt_histo.get_mixture();

    if (cp_ret != NULL)
      {
        ret = new MultivariateMixture(*cp_ret);
        delete cp_ret;
      }
    else
      stat_tool::wrap_util::throw_error(
          "No mixture model available for Mixture Data");
    return ret;
  }

  WRAP_METHOD2(MultivariateMixtureData, extract, DiscreteDistributionData, int, int);
  WRAP_METHOD1(MultivariateMixtureData, extract_marginal, DiscreteDistributionData, int);
  WRAP_METHOD_FILE_ASCII_WRITE( MultivariateMixtureData);
  WRAP_METHOD_PLOT_WRITE( MultivariateMixtureData);
  WRAP_METHOD_SPREADSHEET_WRITE( MultivariateMixtureData);

};

void class_multivariate_mixture_data()
{
  class_< MultivariateMixtureData, bases< Vectors > >
    ("_MultivariateMixtureData", "Multivariate Mixture Data")

    .def("__init__", make_constructor(WRAP::multivariate_mixture_data_from_mixture_data), "Build from a _MultivariateMixture_data object")
    .def(self_ns::str(self))
    //.def("__len__", &MultivariateMixtureData::get_nb_component)
    .def("get_nb_component", &MultivariateMixtureData::get_nb_component, "Return the number of components.")
    DEF_RETURN_VALUE("extract_component", WRAP::extract, args("index"), "Get a particular component for a particular variable. First index is 1")
    DEF_RETURN_VALUE_NO_ARGS("extract_marginal", WRAP::extract_marginal, "Return a _MultivariateMixtureData for a particular variable.")
    DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight, "Return a _MultivariateMixtureData for mixture weights.")
    DEF_RETURN_VALUE_NO_ARGS("extract_mixture", WRAP::extract_mixture, "Return a _MultivariateMixtureData for mixture model")
    .def("file_ascii_write", WRAP::file_ascii_write, "Save _MultivariateMixtureData into a file")
    .def("file_spreadsheet_write", WRAP::spreadsheet_write, "Save _MultivariateMixtureData into a file")
    .def("plot_write", WRAP::plot_write, args("prefix", "title"), "Write GNUPLOT files")
    .def("spreadsheet_write", WRAP::spreadsheet_write, "save data in spreadsheet format")

    ;

/*
    MultivariateMixtureData(const Vectors &vec, int inb_component);
    MultivariateMixtureData(const MultivariateMixture &mixt);
    plotable::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    FrequencyDistribution* get_component(int variable, int index) const { return component[variable][index]; }
*/

};
#undef WRAP


