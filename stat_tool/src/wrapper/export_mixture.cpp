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
#include "stat_tool/markovian.h"
#include "stat_tool/mv_mixture.h"
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


    stat_tool::wrap_util::auto_ptr_array<double>
      weight(new double[nb_component]);

    stat_tool::wrap_util::auto_ptr_array<const Parametric *>
      component(new const Parametric*[nb_component]);

    for(int i=0; i<nb_component; i++)
      {
    weight[i] = boost::python::extract< double >(weights[i]);
    component[i] = boost::python::extract< Parametric *>(dists[i]);
      }

    mix = mixture_building(error, nb_component, weight.get(), component.get());

    if(!mix)
      stat_tool::wrap_util::throw_error(error);


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
      stat_tool::wrap_util::throw_error("Input list cannot be empty");


    stat_tool::wrap_util::auto_ptr_array<const Parametric *>
      component(new const Parametric*[nb_component]);


    for(int i=0; i<nb_component; i++)
      component[i] = boost::python::extract< Parametric *>(dists[i]);

    mix = new Mixture(nb_component, component.get());

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
  static void file_ascii_write(const Mixture m, const char* path, bool exhaustive)
    {
    	bool result = true;
    	Format_error error;

    	result = m.ascii_write(error, path, exhaustive);
    	if (!result)
    	   stat_tool::wrap_util::throw_error(error);
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

    .def("__len__", &Mixture::get_nb_component,
    "Return the number of components") // __len__

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

    .def("file_ascii_write", MixtureWrap::file_ascii_write,
    	     "Save Compound into a file")
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

  static void file_ascii_write(const Mixture_data m, const char* path, bool exhaustive)
  {
  	bool result = true;
  	Format_error error;

  	result = m.ascii_write(error, path, exhaustive);
  	if (!result)
  	   stat_tool::wrap_util::throw_error(error);
    }


};

void class_mixture_data()
{
  class_< Mixture_data, bases< Histogram, STAT_interface > >
    ("_MixtureData",  "Mixture Data")
    .def(self_ns::str(self))
    .def("nb_component", &Mixture_data::get_nb_component,
     "Return the number of components.")
    .def("extract_component", MixtureDataWrap::extract,
     return_value_policy< manage_new_object >(),
     python::arg("index"),
     "Get a particular component. First index is 1")
    .def("extract_weight", MixtureDataWrap::extract_weight,
     return_value_policy< manage_new_object >(),
     "Return a _DistributionData")
    .def("extract_mixture", MixtureDataWrap::extract_mixture,
     return_value_policy< manage_new_object >(),
     "Return a _DistributionData")
     .def("file_ascii_write", MixtureDataWrap::file_ascii_write,
     "Save Compound into a file")
    ;
}



////////////////////////// Class Mv_Mixture //////////////////////////////////


class MvMixtureWrap
{

public:

  static boost::shared_ptr<Mv_Mixture> mv_mixture_from_file(char* filename)
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

  static boost::shared_ptr<Mv_Mixture> mv_mixture_from_mixture(const Mv_Mixture& mixt)
  {
    Mv_Mixture *mix_cp = NULL;

    mix_cp = new Mv_Mixture(mixt, true);

    return boost::shared_ptr<Mv_Mixture>(mix_cp);
  }

  static boost::shared_ptr<Mv_Mixture> mv_mixture_from_components(boost::python::list& weights,
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
    if(nb_component != boost::python::len(dists))
      {
    stat_tool::wrap_util::throw_error("Input lists must have the same length");
      }
    // Test list length
    if(nb_component == 0)
      {
    stat_tool::wrap_util::throw_error("Input lists cannot be empty");
      }

    nb_variable = boost::python::len(dists[0]);

    if(nb_variable == 0)
      {
    stat_tool::wrap_util::throw_error("Number of variable must be positive");
      }


    stat_tool::wrap_util::auto_ptr_array<double>
      weight(new double[nb_component]);

    stat_tool::wrap_util::auto_ptr_array<Parametric_process *>
      pcomponent(new Parametric_process*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<Nonparametric_process *>
      npcomponent(new Nonparametric_process*[nb_variable]);

    stat_tool::wrap_util::auto_ptr_array<bool>
      is_parametric(new bool[nb_variable]);

    int i, var;

    for(i=0; i<nb_component; i++)
      weight[i] = boost::python::extract< double >(weights[i]);

    for (var=0; var < nb_variable; var++) {

      Parametric **pprocess = new Parametric*[nb_component];
      Distribution **npprocess = new Distribution*[nb_component];

      pprocess[0] = NULL;
      npprocess[0] = NULL;

      for(i=0; i<nb_component; i++) {
    boost::python::extract< Parametric* > x(dists[i][var]);
    if (x.check()) {
      pprocess[i] = x();
      npprocess[i] = NULL;
      if (i == 0)
        is_parametric[var] = true;
      else {
        if (!is_parametric[var])
          stat_tool::wrap_util::throw_error("Distributions must be the same type for a given variable");
      }
    }
    else {
      boost::python::extract< Distribution *> x(dists[i][var]);
      if (x.check()) {
        npprocess[i] = x();
        pprocess[i] = NULL;
        if (i == 0)
          is_parametric[var] = false;
        else {
          if (is_parametric[var])
        stat_tool::wrap_util::throw_error("Distributions must be the same type for a given variable");
        }
      }
      else
        stat_tool::wrap_util::throw_error("Bad type for list element: must be Parametric or Distribution");

    }
      } // end for (i)
      if (pprocess[0] != NULL) {
    pcomponent[var] = new Parametric_process(nb_component, pprocess);
    npcomponent[var] = NULL;
      }
      else {
    pcomponent[var] = NULL;
    npcomponent[var] = new Nonparametric_process(nb_component, npprocess);
      }

      /*for(i=0; i<nb_component; i++) {
    if (npprocess[i] != NULL)
      delete npprocess[i];
    if (pprocess[i] != NULL)
    delete pprocess[i];
    }*/

      delete [] pprocess;
      delete [] npprocess;
    } // end for (var)

    mix = mv_mixture_building(error, nb_component, nb_variable,
                  weight.get(), pcomponent.get(), npcomponent.get());

    for(var=0; var<nb_variable; var++) {
      if (pcomponent[var] != NULL)
      delete pcomponent[var];
      if (npcomponent[var] != NULL)
      delete npcomponent[var];
    }

    if(mix == NULL)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<Mv_Mixture>(mix);
  }


  static Mv_Mixture_data* simulation(const Mv_Mixture& mixt, int nb_element)
  {
    Format_error error;
    Mv_Mixture_data* ret = NULL;

    ret = mixt.simulation(error, nb_element);
    if(ret == NULL) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Parametric_model* extract_weight(const Mv_Mixture& mixt)
  {
    Parametric_model* ret;
    Mv_Mixture_data* mixt_data = NULL;

    mixt_data = mixt.get_mixture_data();
    ret = new Parametric_model(*(mixt.get_weight()),
                   (mixt_data ? mixt_data->get_weight() : NULL));
    return ret;
  }


  static Parametric_model* extract_mixture(const Mv_Mixture& mixt, int ivariable)
  {
    Format_error error;
    Parametric_model* ret;
    Mv_Mixture_data* mixt_data = NULL;
    Distribution *marginal = mixt.extract_distribution(error, ivariable);
    Histogram *marginal_hist = NULL;

    if (marginal == NULL) stat_tool::wrap_util::throw_error(error);

    mixt_data = mixt.get_mixture_data();

    if (mixt_data != NULL)
      marginal_hist = mixt_data->get_marginal(ivariable);

    ret = new Parametric_model(*marginal, marginal_hist);

    if (marginal_hist != NULL)
      delete marginal_hist;

    if (marginal != NULL) {
      delete marginal;
      marginal = NULL;
    }

    return ret;
  }

  static Mv_Mixture_data* extract_data(const Mv_Mixture& mixt)
  {
    Format_error error;
    Mv_Mixture_data* ret = NULL;

    ret = mixt.extract_data(error);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static bool _is_parametric(const Mv_Mixture& mixt, int ivariable)
  {
    ostringstream error_message;

    if ((ivariable < 0) || (ivariable > mixt.get_nb_component())) {
      error_message << "Bad variable index: " << ivariable;
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
    }
    else
      return mixt.is_parametric(ivariable);
  }



  static void file_ascii_write(const Mv_Mixture& m, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = m.ascii_write(error, path, exhaustive);
    if (!result)
       stat_tool::wrap_util::throw_error(error);

  }

  static void plot_write(const Mv_Mixture& m, const char* path,  const char* title)
  {
    bool result = true;
    Format_error error;

    result = m.plot_write(error, path, title);
    if (!result)
      stat_tool::wrap_util::throw_error(error);
  }

  static void state_permutation(const Mv_Mixture& mix,
                                boost::python::list perm)
   {
      bool status = true, several_errors = false;
      int llength, i;
      ostringstream error_message;
      object o;
      Format_error error;
      int *iperm = NULL;

      llength= boost::python::len(perm);
      if ((llength > 0) && (llength == mix.get_nb_component()))
      {
         iperm = new int[llength];
         for (i= 0; i < llength; i++)
         {
            o= perm[i];
            try
            {
               extract<int> x(o);
               if (x.check())
                  iperm[i]= x();
               else
                  status= false;
            }
            catch (...)
            {
               status= false;
            }
            if (!status)
            {
               if (several_errors)
                  error_message << endl;
               else
                  several_errors= true;
               error_message << "incorrect type for element " << i
                             << " of argument list: expecting an int ";
            }
         }
         if (!status)
         {
            delete [] iperm;
            iperm= NULL;
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
         delete [] iperm;
         iperm= NULL;
      }
   }

};



// Boost declaration

void class_mv_mixture()
{
  class_< Mv_Mixture, bases< STAT_interface > >
    ("_MvMixture", "Multivariate Mixture Distribution")

    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_file),
     "Build from a filename"
     )

    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_mixture),
     "Build from a _MvMixture object"
     )

    .def("__init__", make_constructor(MvMixtureWrap::mv_mixture_from_components),
     "Build from a list of weights and a list of list of distributions\n"
     "(components and variables)")

    .def(self_ns::str(self)) // __str__
    .def("__len__", &Mv_Mixture::get_nb_component,
    "Return the number of components") // __len__

    .def("simulate", MvMixtureWrap::simulation,
     return_value_policy< manage_new_object >(),
     python::arg("nb_element"),
     "Simulate nb_element elements")

    .def("nb_component", &Mv_Mixture::get_nb_component,
     "Return the number of components")

    .def("nb_variable", &Mv_Mixture::get_nb_variable,
     "Return the number of variables")

    .def("extract_weight", MvMixtureWrap::extract_weight,
     return_value_policy< manage_new_object >(),
     "Return the weight distribution")

    .def("extract_mixture", MvMixtureWrap::extract_mixture,
     return_value_policy< manage_new_object >(),
     "Return the _MvMixture distribution")

    .def("extract_data", MvMixtureWrap::extract_data,
     return_value_policy< manage_new_object >(),
     "Return the associated _MvMixtureData object"
     )

    .def("_is_parametric", MvMixtureWrap::_is_parametric,
     python::args("variable"),
     "Return True if the variable is parametric"
     )

    .def("file_ascii_write", MvMixtureWrap::file_ascii_write,
     "Save _MvMixture into a file")

    .def("plot_write", MvMixtureWrap::plot_write,
     python::args("prefix", "title"),
     "Write GNUPLOT files")

    .def("state_permutation", MvMixtureWrap::state_permutation,
     "permutation of the model states")

     ;

}



////////////////////////// Class Mv_Mixture_data //////////////////////////////////

class MvMixtureDataWrap
{

public:

  static boost::shared_ptr<Mv_Mixture_data> mv_mixture_data_from_mixture_data(const Mv_Mixture_data& mixt)
  {
    Mv_Mixture_data *mix_cp = NULL;

    mix_cp = new Mv_Mixture_data(mixt, true);

    return boost::shared_ptr<Mv_Mixture_data>(mix_cp);
  }

  static Distribution_data* extract(const Mv_Mixture_data& d, int ivariable, int index)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = d.extract(error, ivariable, index);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data* extract_marginal(const Mv_Mixture_data& d, int ivariable)
  {
    Format_error error;
    Distribution_data* ret = NULL;

    ret = d.extract_marginal(error, ivariable);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data* extract_weight(const Mv_Mixture_data& mixt_histo)
  {
    Distribution_data* ret;

    ret = new Distribution_data(*(mixt_histo.get_weight()) ,
                mixt_histo.get_mixture()->get_weight());

    return ret;
  }


  static Mv_Mixture* extract_mixture(const Mv_Mixture_data& mixt_histo)
  {
    Mv_Mixture* ret, *cp_ret = NULL;

    cp_ret = mixt_histo.get_mixture();

    if (cp_ret != NULL) {
      ret = new Mv_Mixture(*cp_ret);
      delete cp_ret;
    }
    else
      stat_tool::wrap_util::throw_error("No mixture model available for Mixture Data");
    return ret;
  }

  static void file_ascii_write(const Mv_Mixture_data& m, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = m.ascii_write(error, path, exhaustive);
    if (!result)
       stat_tool::wrap_util::throw_error(error);

  }

  static void plot_write(const Mv_Mixture_data& m, const char* path,
             const char* title)
  {
    bool result = true;
    Format_error error;

    result = m.plot_write(error, path, title);
    if (!result)
      stat_tool::wrap_util::throw_error(error);
  }

};

void class_mv_mixture_data()
{
  class_< Mv_Mixture_data, bases< Vectors > >
    ("_MvMixtureData", "Multivariate Mixture Data")

    .def("__init__", make_constructor(MvMixtureDataWrap::mv_mixture_data_from_mixture_data),
     "Build from a _MvMixture_data object"
     )

    .def(self_ns::str(self))

    .def("nb_component", &Mv_Mixture_data::get_nb_component,
     "Return the number of components."
     )

    .def("extract_component", MvMixtureDataWrap::extract,
     return_value_policy< manage_new_object >(),
     python::arg("index"),
     "Get a particular component for a particular variable. First index is 1"
     )

    .def("extract_marginal", MvMixtureDataWrap::extract_marginal,
     return_value_policy< manage_new_object >(),
     "Return a _MvMixtureData for a particular variable."
     )

    .def("extract_weight", MvMixtureDataWrap::extract_weight,
     return_value_policy< manage_new_object >(),
     "Return a _MvMixtureData for mixture weights."
     )

    .def("extract_mixture", MvMixtureDataWrap::extract_mixture,
     return_value_policy< manage_new_object >(),
     "Return a _MvMixtureData for mixture model"
     )

    .def("file_ascii_write", MvMixtureDataWrap::file_ascii_write,
     "Save _MvMixtureData into a file")

    .def("plot_write", MvMixtureDataWrap::plot_write,
     python::args("prefix", "title"),
     "Write GNUPLOT files")
    ;
}
