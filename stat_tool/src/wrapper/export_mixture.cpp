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

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/convolution.h"
#include "stat_tool/mixture.h"
#include "stat_tool/markovian.h"
#include "stat_tool/mv_mixture.h"
#include "stat_tool/compound.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;



////////////////////// Export Mixture ////////////////////////////////////////
#define WRAP MixtureWrap
class WRAP
{

public:

  static boost::shared_ptr<Mixture>
  mixture_from_file(char* filename)
  {
    Format_error error;
    Mixture *mix = NULL;
    mix = mixture_ascii_read(error, filename);
    if (!mix)
      stat_tool::wrap_util::throw_error(error);
    return boost::shared_ptr<Mixture>(mix);
  }

  static boost::shared_ptr<Mixture>
  mixture_from_dists(boost::python::list& weights, boost::python::list& dists)
  {
    Format_error error;
    Mixture *mix = NULL;
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

    stat_tool::wrap_util::auto_ptr_array<const Parametric*> component(
        new const Parametric*[nb_component]);

    int i = 0;
    for (i = 0; i < nb_component; i++)
      {
        weight[i] = boost::python::extract<double>(weights[i]);

        boost::python::extract<Parametric*> get_param(dists[i]);
        if (get_param.check())
          {
            component[i] = new Parametric(*get_param());
          }
        else
          {
            component[i] = new Parametric(
                *boost::python::extract<Distribution*>(dists[i])());
          }

      }

    mix = mixture_building(error, nb_component, weight.get(), component.get());

    if (!mix)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<Mixture>(mix);
  }

  static boost::shared_ptr<Mixture>
  mixture_from_unknown_component(boost::python::list& dists)
  {
    Format_error error;
    Mixture *mix = NULL;
    int nb_component = 0;

    nb_component = boost::python::len(dists);

    // Test list length
    if (nb_component == 0)
      stat_tool::wrap_util::throw_error("Input list cannot be empty");

    stat_tool::wrap_util::auto_ptr_array<const Parametric *> component(
        new const Parametric*[nb_component]);

    for (int i = 0; i < nb_component; i++)
      {
        component[i] = boost::python::extract<Parametric *>(dists[i]);
      }

    mix = new Mixture(nb_component, component.get());
    if (!mix)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<Mixture>(mix);
  }

  static Parametric_model*
  extract_weight(const Mixture& mixt)
  {
    Format_error error;
    Parametric_model* ret;
    Mixture_data* mixt_data = NULL;

    mixt_data = mixt.get_mixture_data();
    ret = new Parametric_model(*(mixt.get_weight()),
        (mixt_data ? mixt_data->get_weight() : NULL));
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static Parametric_model*
  extract_mixture(Mixture& mixture_input)
  {
    Format_error error;
    Parametric_model* ret;
    Mixture_data* mixture_data = NULL;

    mixture_data = mixture_input.get_mixture_data();

    ret = new Parametric_model(*((Distribution*) (&mixture_input)),
        (Histogram*) mixture_data);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Parametric_model*
  extract_component(const Mixture& input, int var1)
  {
    Format_error error;
    Parametric_model* ret = NULL;

    ret = input.extract(error, var1);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // component case
  static Mixture_data*
  extract_data(const Mixture& input)
  {
    Format_error error;
    Mixture_data* ret = NULL;

    ret = input.extract_data(error);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  WRAP_METHOD1(Mixture, simulation, Mixture_data, int);
  WRAP_METHOD_FILE_ASCII_WRITE( Mixture);
  WRAP_METHOD_SPREADSHEET_WRITE( Mixture);

  static MultiPlotSet*
  get_plotable(const Mixture& mixt)
  {
    Format_error error;
    MultiPlotSet* ret = mixt.get_plotable();
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
    return ret;
  }

  static MultiPlotSet*
  survival_get_plotable(const Mixture& p)
  {
    Format_error error;
    MultiPlotSet* ret = p.survival_get_plotable(error);
    if (!ret)
      ERROR;
    return ret;
  }

};



// Boost declaration
void class_mixture()
{

  class_< Mixture, bases< Distribution, STAT_interface > >
  ("_Mixture", "Mixture Distribution")
  // constructors

  .def("__init__", make_constructor(MixtureWrap::mixture_from_file),
      "Build from a filename" )
  .def("__init__", make_constructor(MixtureWrap::mixture_from_dists),
      "Build from a list of weights and a list of distributions")
  .def("__init__", make_constructor(MixtureWrap::mixture_from_unknown_component),
      "Build from unknown components") // internal use

  // Python Operators
  .def("__len__", &Mixture::get_nb_component,
      "Return the number of components") // __len__
  .def(self_ns::str(self)) // __str__

  // properties
  .add_property("nb_component", &Mixture::get_nb_component,
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
      "Return the Mixture distribution")
  DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,
      "Return the associated _MixtureData object")


  //others
  //  DEF_RETURN_VALUE_NO_ARGS("get_plotable", &STAT_interface::get_plotable,"Return a plotable (no parameters)");
  .def("file_ascii_write", WRAP::file_ascii_write,
      "Save Compound into a file")
  .def("spreadsheet_write", WRAP::spreadsheet_write,
      "save data in spreadsheet format")
  DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable,
      "return plotable")
  DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable,
      "Return a survival plotable")

  ;

  /*

  Mixture();
  Mixture(const Mixture &mixt , bool *component_flag , int inb_value);
  Mixture(int inb_component , const Parametric **pcomponent);
  Mixture(const Mixture &mixt , bool data_flag = true)    :Distribution(mixt)
    { copy(mixt , data_flag); }


  void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD ,
                     bool component_flag = true);
  double likelihood_computation(const Mixture_data &mixt_histo) const;

  // acces membres de la classe

  Parametric* get_component(int index) const { return component[index]; }
*/
}
#undef WRAP


////////////////////////// Class Mixture_data //////////////////////////////////
#define WRAP MixtureDataWrap
class MixtureDataWrap
{

public:

  WRAP_METHOD1(Mixture_data, extract, Distribution_data, int);
  WRAP_METHOD_FILE_ASCII_WRITE( Mixture_data);

  static Distribution_data*
  extract_weight(const Mixture_data& mixt_histo)
  {
    Distribution_data* ret;
    ret = new Distribution_data(*(mixt_histo.get_weight()),
        mixt_histo.get_mixture()->get_weight());
    return ret;
  }

  static Distribution_data*
  extract_mixture(const Mixture_data& mixt_histo)
  {
    Distribution_data* ret;
    ret = new Distribution_data(mixt_histo, mixt_histo.get_mixture());
    return ret;
  }

};

void class_mixture_data()
{
  class_< Mixture_data, bases< Histogram, STAT_interface > >
  ("_MixtureData",  "Mixture Data")

  // Python Operators
  .def(self_ns::str(self)) //str

  // properties
 .add_property("nb_component", &Mixture_data::get_nb_component,
      "Return the number of components.")

  // getters
  DEF_RETURN_VALUE("get_component", &Mixture_data::get_component,
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
    Mixture_data(const Histogram &histo , int inb_component);
    Mixture_data(const Mixture &mixt);
    Mixture_data(const Mixture_data &mixt_histo , bool model_flag = true) :Histogram(mixt_histo) { copy(mixt_histo , model_flag); }

    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix , const char *title = 0) const;
    plotable::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    Histogram* get_component(int index) const { return component[index]; }
    */
    ;
}

#undef WRAP



////////////////////////// Class Mv_Mixture //////////////////////////////////

#define WRAP MvMixtureWrap
class MvMixtureWrap
{

public:

  static boost::shared_ptr<Mv_Mixture>
  mv_mixture_from_file(char* filename)
  {
    Format_error error;
    Mv_Mixture *mix = NULL;
    mix = mv_mixture_ascii_read(error, filename);

    if (mix == NULL)
      {
        stat_tool::wrap_util::throw_error(error);
      }

    return boost::shared_ptr<Mv_Mixture>(mix);
  }

  static boost::shared_ptr<Mv_Mixture>
  mv_mixture_from_mixture(const Mv_Mixture& mixt)
  {
    Mv_Mixture *mix_cp = NULL;

    mix_cp = new Mv_Mixture(mixt, true);

    return boost::shared_ptr<Mv_Mixture>(mix_cp);
  }

  static boost::shared_ptr<Mv_Mixture>
  mv_mixture_from_components(boost::python::list& weights,
      boost::python::list& dists)
  {
    Format_error error;
    Mv_Mixture *mix = NULL;
    int nb_component = 0;
    int nb_variable = 0;
    // Parametric *pcomp = NULL;
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

    stat_tool::wrap_util::auto_ptr_array<Parametric_process *> pcomponent(
        new Parametric_process*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<Nonparametric_process *> npcomponent(
        new Nonparametric_process*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<bool> is_parametric(
        new bool[nb_variable]);

    int i, var;

    for (i = 0; i < nb_component; i++)
      weight[i] = boost::python::extract<double>(weights[i]);

    for (var = 0; var < nb_variable; var++)
      {

        Parametric **pprocess = new Parametric*[nb_component];
        Distribution **npprocess = new Distribution*[nb_component];

        pprocess[0] = NULL;
        npprocess[0] = NULL;

        for (i = 0; i < nb_component; i++)
          {
            boost::python::extract<Parametric*> x(dists[i][var]);
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
                      "Bad type for list element: must be Parametric or Distribution");

              }
          } // end for (i)
        if (pprocess[0] != NULL)
          {
            pcomponent[var] = new Parametric_process(nb_component, pprocess);
            npcomponent[var] = NULL;
          }
        else
          {
            pcomponent[var] = NULL;
            npcomponent[var] = new Nonparametric_process(nb_component,
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

    mix = mv_mixture_building(error, nb_component, nb_variable, weight.get(),
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

    return boost::shared_ptr<Mv_Mixture>(mix);
  }

  static Parametric_model*
  extract_weight(const Mv_Mixture& mixt)
  {
    Parametric_model* ret;
    Mv_Mixture_data* mixt_data = NULL;

    mixt_data = mixt.get_mixture_data();
    ret = new Parametric_model(*(mixt.get_weight()),
        (mixt_data ? mixt_data->get_weight() : NULL));
    return ret;
  }

  static Parametric_model*
  extract_mixture(const Mv_Mixture& mixt, int ivariable)
  {
    Format_error error;
    Parametric_model* ret;
    Mv_Mixture_data* mixt_data = NULL;
    Distribution *marginal = mixt.extract_distribution(error, ivariable);
    Histogram *marginal_hist = NULL;

    if (marginal == NULL)
      stat_tool::wrap_util::throw_error(error);

    mixt_data = mixt.get_mixture_data();

    if (mixt_data != NULL)
      marginal_hist = mixt_data->get_marginal(ivariable);

    ret = new Parametric_model(*marginal, marginal_hist);

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

  static bool
  _is_parametric(const Mv_Mixture& mixt, int ivariable)
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
  state_permutation(const Mv_Mixture& mix, boost::python::list perm)
  {
    bool status = true, several_errors = false;
    int llength, i;
    ostringstream error_message;
    object o;
    Format_error error;
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

    WRAP_METHOD0(Mv_Mixture, extract_data, Mv_Mixture_data);
    WRAP_METHOD_FILE_ASCII_WRITE( Mv_Mixture);
    WRAP_METHOD_PLOT_WRITE( Mv_Mixture);
    WRAP_METHOD_SPREADSHEET_WRITE( Mv_Mixture);
    WRAP_METHOD1(Mv_Mixture, simulation, Mv_Mixture_data, int);

};



// Boost declaration



void class_mv_mixture()
{
  class_< Mv_Mixture, bases< STAT_interface > >
    ("_MvMixture", "Multivariate Mixture Distribution")

    // constructors
    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_file),
        "Build from a filename")
    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_mixture),
        "Build from a _MvMixture object")
    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_components),
        "Build from a list of weights and a list of list of distributions\n" "(components and variables)")

    // Python Operators
    .def(self_ns::str(self)) // __str__
    .def("__len__", &Mv_Mixture::get_nb_component,
        "Return the number of components") // __len__

    // properties
    .add_property("nb_variable", &Mv_Mixture::get_nb_variable,
        "Return the number of variables")
    .add_property("nb_component", &Mv_Mixture::get_nb_component,
        "Return the number of variables")

    // return object and no arguments needed
    DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight,
        "Return the weight distribution")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data,
        "Return the associated _MvMixtureData object")


    // return object and arguments needed
    DEF_RETURN_VALUE("simulate", WRAP::simulation,
        args("nb_element"), "simulate(self, nb_element) -> _MvMixtureData. \n\n" "Simulate nb_element elements")
    DEF_RETURN_VALUE("extract_mixture", WRAP::extract_mixture,
        args("variable"), "extract_mixture(self, variable) -> _Distribution. \n\n" "Return the _MvMixture distribution")


    // no object returned, args required
    .def("_is_parametric", WRAP::_is_parametric,
        args("variable"),"_is_parametric(self, variable) -> bool. \n\n" "Return True if the variable is parametric")
    .def("file_ascii_write", WRAP::file_ascii_write,
        args("path", "exhaustive_flag"), "file_ascii_write(self, path, exhaustive_flag) -> None. \n\n""Save _MvMixture into a file")
    .def("plot_write", WRAP::plot_write,
        args("prefix", "title"),"plot_write(self, prefix, title) -> None. \n\n" "Write GNUPLOT files")

    // no object returned, no arguments required
    .def("state_permutation", WRAP::state_permutation,
        "permutation of the model states")
    .def("spreadsheet_write", WRAP::spreadsheet_write,
        "save data in spreadsheet format")

     ;

    /*

      Mv_Mixture(int inb_component , double *pweight , int inb_variable, Parametric_process **ppcomponent, Nonparametric_process **pnpcomponent);
      Mv_Mixture(int inb_component , int inb_variable, const Parametric_process **ppcomponent,  const Nonparametric_process **pnpcomponent);
      Mv_Mixture(const Mv_Mixture &mixt , bool *variable_flag , int inb_variable);
      Mv_Mixture(int inb_component, int inb_variable, int *nb_value, bool *force_param=NULL);

      Parametric_model* extract_parametric_model(Format_error &error , int ivariable,
                             int index) const;
      Distribution* extract_nonparametric_model(Format_error &error , int ivariable,
                            int index) const;

      plotable::MultiPlotSet* get_plotable() const;

      double likelihood_computation(const Vectors &mixt_data, bool log_computation=false) const;

      Mv_Mixture_data* cluster(Format_error &error,  const Vectors &vec,int algorithm=VITERBI) const;

      Mv_Mixture_data* get_mixture_data() const { return mixture_data; }
      Parametric_process* get_parametric_process(int variable) const;
      Nonparametric_process* get_nonparametric_process(int variable) const;
      Parametric* get_parametric_component(int variable, int index) const;
      Distribution* get_nonparametric_component(int variable, int index) const;*/
}

#undef WRAP

////////////////////////// Class Mv_Mixture_data //////////////////////////////////

#define WRAP MvMixtureDataWrap
class WRAP
{

public:

  static boost::shared_ptr<Mv_Mixture_data>
  mv_mixture_data_from_mixture_data(const Mv_Mixture_data& mixt)
  {
    Mv_Mixture_data *mix_cp = NULL;

    mix_cp = new Mv_Mixture_data(mixt, true);

    return boost::shared_ptr<Mv_Mixture_data>(mix_cp);
  }

  static Distribution_data*
  extract_weight(const Mv_Mixture_data& mixt_histo)
  {
    Distribution_data* ret;

    ret = new Distribution_data(*(mixt_histo.get_weight()),
        mixt_histo.get_mixture()->get_weight());

    return ret;
  }

  static Mv_Mixture*
  extract_mixture(const Mv_Mixture_data& mixt_histo)
  {
    Mv_Mixture* ret, *cp_ret = NULL;

    cp_ret = mixt_histo.get_mixture();

    if (cp_ret != NULL)
      {
        ret = new Mv_Mixture(*cp_ret);
        delete cp_ret;
      }
    else
      stat_tool::wrap_util::throw_error(
          "No mixture model available for Mixture Data");
    return ret;
  }
  WRAP_METHOD2(Mv_Mixture_data, extract, Distribution_data, int, int);
  WRAP_METHOD1(Mv_Mixture_data, extract_marginal, Distribution_data, int);
  WRAP_METHOD_FILE_ASCII_WRITE( Mv_Mixture_data);
  WRAP_METHOD_PLOT_WRITE( Mv_Mixture_data);
  WRAP_METHOD_SPREADSHEET_WRITE( Mv_Mixture_data);

};

void class_mv_mixture_data()
{
  class_< Mv_Mixture_data, bases< Vectors > >
    ("_MvMixtureData", "Multivariate Mixture Data")

    .def("__init__", make_constructor(WRAP::mv_mixture_data_from_mixture_data), "Build from a _MvMixture_data object")
    .def(self_ns::str(self))
    //.def("__len__", &Mv_Mixture_data::get_nb_component)
    .def("get_nb_component", &Mv_Mixture_data::get_nb_component, "Return the number of components.")
    DEF_RETURN_VALUE("extract_component", WRAP::extract, args("index"), "Get a particular component for a particular variable. First index is 1")
    DEF_RETURN_VALUE_NO_ARGS("extract_marginal", WRAP::extract_marginal, "Return a _MvMixtureData for a particular variable.")
    DEF_RETURN_VALUE_NO_ARGS("extract_weight", WRAP::extract_weight, "Return a _MvMixtureData for mixture weights.")
    DEF_RETURN_VALUE_NO_ARGS("extract_mixture", WRAP::extract_mixture, "Return a _MvMixtureData for mixture model")
    .def("file_ascii_write", WRAP::file_ascii_write, "Save _MvMixtureData into a file")
    .def("file_spreadsheet_write", WRAP::spreadsheet_write, "Save _MvMixtureData into a file")
    .def("plot_write", WRAP::plot_write, args("prefix", "title"), "Write GNUPLOT files")
    .def("spreadsheet_write", WRAP::spreadsheet_write, "save data in spreadsheet format")

    ;

/*
    Mv_Mixture_data(const Vectors &vec, int inb_component);
    Mv_Mixture_data(const Mv_Mixture &mixt);
    plotable::MultiPlotSet* get_plotable() const;

    double information_computation() const;

    Histogram* get_component(int variable, int index) const { return component[variable][index]; }
*/

};
#undef WRAP


