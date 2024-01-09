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
 *        $Id: export_frequency_distribution.cpp 18035 2015-04-23 07:15:42Z guedon $
 *
 *-----------------------------------------------------------------------------*/



#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include "stat_tool/distribution.h"
#include "stat_tool/convolution.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/compound.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;



////////////////////// Export FrequencyDistribution ////////////////////////////////////////


// Wrapper class class
class FrequencyDistributionWrap
{

public:

  // __len__
  static int
  histo_get_nb_value(FrequencyDistribution &h)
  {
    return h.nb_value;
  }

  static int
  histo_get_nb_element(FrequencyDistribution &h)
  {
    return h.nb_element;
  }

  // __getitem__
  static int
  histo_get_item(FrequencyDistribution &h, const object& key)
  {
    int val = extract<int> (key);

    if ((val >= 0) && (val < h.nb_value))
      return h.frequency[val];
    else
      {
        PyErr_SetString(PyExc_IndexError, "Bad Index");
        throw_error_already_set();
      }
  }

  static std::string
  f_comparison(const FrequencyDistribution &h, const FrequencyDistribution &histo)
  {
    ostringstream output;
    h.F_comparison(output, histo);
    return output.str();

  }

  static std::string
  t_comparison(const FrequencyDistribution &h, const FrequencyDistribution &histo)
  {
    ostringstream output;
    h.t_comparison(output, histo);
    return output.str();
  }

  static std::string
  wmw_comparison(const FrequencyDistribution &h, const FrequencyDistribution &histo)
  {
    ostringstream output;
    StatError error;
    if (!h.wilcoxon_mann_whitney_comparison(error, output, histo))
      {
        stat_tool::wrap_util::throw_error(error);
      }
    return output.str();

  }

  static std::string
  comparison(const FrequencyDistribution &h, boost::python::tuple& histos, int type,
      const char* filename, char format)
  {
    ostringstream output;
    StatError error;
    int nb_histo = boost::python::len(histos);

    // Test list length
    if (nb_histo == 0)
      stat_tool::wrap_util::throw_error("No frequency distribution to compare");

    stat_tool::wrap_util::auto_ptr_array<const FrequencyDistribution*> ihisto(
        new const FrequencyDistribution*[nb_histo]);

    // Extract list element
    for (int i = 0; i < nb_histo; i++)
      ihisto[i] = boost::python::extract<FrequencyDistribution*>(histos[i]);

    // call comparaison
    bool res = h.comparison(error, output, nb_histo, ihisto.get(), type,
        filename, format);

    if (!res)
      stat_tool::wrap_util::throw_error(error);

    return output.str();
  }

  // Estimation functions

  static DiscreteMixture*
  discrete_mixture_estimation1(const FrequencyDistribution &h, const DiscreteMixture& imixt,
      boost::python::list& estimate_tuple, int min_inf_bound, bool flag,
      bool component_flag)
  {
    DiscreteMixture* ret;
    StatError error;
    bool estimate[DISCRETE_MIXTURE_NB_COMPONENT];

    //Parse estimate list
    int nb_element = boost::python::len(estimate_tuple);

    for (int i = 0; i < nb_element; i++)
      estimate[i] = boost::python::extract<bool>(estimate_tuple[i]);

    ret = h.discrete_mixture_estimation(error, imixt, estimate, min_inf_bound, flag,
        component_flag);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;

  }

  static DiscreteMixture*
  discrete_mixture_estimation2(const FrequencyDistribution &h, boost::python::list& ident_list,
      int min_inf_bound, bool flag, bool component_flag, int penalty)
  {
    DiscreteMixture* ret;
    ostringstream output;
    StatError error;

    int nb_component = boost::python::len(ident_list);
    int ident[DISCRETE_MIXTURE_NB_COMPONENT];

    for (int i = 0; i < nb_component; i++)
      ident[i] = boost::python::extract<int>(ident_list[i]);

    ret = h.discrete_mixture_estimation(error, output, 1, nb_component, ident,
        min_inf_bound, flag, component_flag, penalty);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    cout << output.str() << endl;

    return ret;
  }

  static Convolution*
  convolution_estimation1(const FrequencyDistribution &h, const DiscreteParametric &known_dist,
      const DiscreteParametric &unknown_dist, int estimator, int nb_iter,
      double weight, int penalty_type, int outside)
  {
    Convolution* ret;
    ostringstream output;
    StatError error;

    ret = h.convolution_estimation(error, output, known_dist, unknown_dist,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Convolution*
  convolution_estimation2(const FrequencyDistribution &h, const DiscreteParametric &known_dist,
      int min_inf_bound, int estimator, int nb_iter, double weight,
      int penalty_type, int outside)
  {
    Convolution* ret;
    ostringstream output;
    StatError error;

    ret = h.convolution_estimation(error, output, known_dist, min_inf_bound,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Compound*
  compound_estimation1(const FrequencyDistribution &h, const DiscreteParametric &sum_dist,
      const DiscreteParametric &dist, char type, int estimator, int nb_iter,
      double weight, int penalty_type, int outside)
  {
    Compound* ret;
    ostringstream output;
    StatError error;

    ret = h.compound_estimation(error, output, sum_dist, dist, type, estimator,
        nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Compound*
  compound_estimation2(const FrequencyDistribution &h, const DiscreteParametric &known_dist,
      char type, int min_inf_bound, int estimator, int nb_iter, double weight,
      int penalty_type, int outside)
  {
    Compound* ret;
    ostringstream output;
    StatError error;

    ret = h.compound_estimation(error, output, known_dist, type, min_inf_bound,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Merge FrequencyDistribution in a distribution_data
  static DiscreteDistributionData*
  merge_frequency_distributions(const FrequencyDistribution &h, const boost::python::list& histo_list)
  {
    DiscreteDistributionData *histo = NULL;
    int nb_element = boost::python::len(histo_list) + 1;

    stat_tool::wrap_util::auto_ptr_array<const FrequencyDistribution*> pelement(
        new const FrequencyDistribution*[nb_element]);

    pelement[0] = &h;

    for (int i = 1; i < nb_element; i++)
      pelement[i] = extract<FrequencyDistribution*> (histo_list[i - 1]);

    // create new distribution_data object
    histo = new DiscreteDistributionData(nb_element, pelement.get());
    if (!histo)
      stat_tool::wrap_util::throw_error(
          "Could not initialize FrequencyDistribution from arguments");

    return histo;
  }

  WRAP_METHOD3(FrequencyDistribution, parametric_estimation, DiscreteParametricModel, int, int, bool);
  WRAP_METHOD3(FrequencyDistribution, value_select, DiscreteDistributionData, int, int, bool);
  WRAP_METHOD1(FrequencyDistribution, shift, DiscreteDistributionData, int);
  WRAP_METHOD1(FrequencyDistribution, fit, DiscreteParametricModel, DiscreteParametric);

  // Cluster
  static DiscreteDistributionData*
  cluster_step(const FrequencyDistribution &h, int step)
  {
    StatError error;
    DiscreteDistributionData* ret = h.cluster(error, step);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static DiscreteDistributionData*
  cluster_information(const FrequencyDistribution &h, float ratio)
  {
    StatError error;
    std::stringstream output;

    DiscreteDistributionData* ret = h.cluster(error, ratio, output);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static DiscreteDistributionData*
  cluster_limit(const FrequencyDistribution &h, boost::python::list& limit)
  {

    StatError error;

    int nb_limit = len(limit);
    stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_limit]);

    for (int i = 0; i < nb_limit; i++)
      l[i] = extract<int> (limit[i]);

    DiscreteDistributionData* ret = h.cluster(error, nb_limit + 1, l.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static DiscreteDistributionData*
  transcode(const FrequencyDistribution &h, boost::python::list& symbol)
  {

    StatError error;

    int nb_symbol = len(symbol);
    stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

    int expected_nb_symbol = h.nb_value - h.offset;
    if (nb_symbol != expected_nb_symbol)
      stat_tool::wrap_util::throw_error("Bad number of Symbol");

    for (int i = 0; i < nb_symbol; i++)
      l[i] = extract<int> (symbol[i]);

    DiscreteDistributionData* ret = h.transcode(error, l.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static MultiPlotSet*
  get_plotable(const FrequencyDistribution &p, const boost::python::list& hist_list)
  {
    StatError error;
    int nb_hist = boost::python::len(hist_list);
    stat_tool::wrap_util::auto_ptr_array<const FrequencyDistribution *> hists(
        new const FrequencyDistribution*[nb_hist]);

    for (int i = 0; i < nb_hist; i++)
      hists[i] = extract<const FrequencyDistribution*> (hist_list[i]);

    const FrequencyDistribution** d = hists.get();

    MultiPlotSet* ret = p.get_plotable_frequency_distributions(error, nb_hist, d);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

};



// Boost declaration

void class_frequency_distribution()
{

   // _FrequencyDistribution
   class_<FrequencyDistribution> ("_FrequencyDistribution", no_init)
  .def(self == self)
  .def(self != self) .def(self_ns::str(self))
  .def("__len__", FrequencyDistributionWrap::histo_get_nb_element)
  .def("__getitem__", FrequencyDistributionWrap::histo_get_item)

  // Comparison
  .def("compare", FrequencyDistributionWrap::comparison,
      args("frequency distribution"),"Comparison of frequency distributions")
  .def("f_comparison", FrequencyDistributionWrap::f_comparison,
      args("frequency distribution"), "F comparison of frequency distributions")
  .def("t_comparison", FrequencyDistributionWrap::t_comparison, args("frequency distribution"),
      "T comparison of frequency distributions")
  .def("wmw_comparison", FrequencyDistributionWrap::wmw_comparison,
      args("frequency distribution"),
      " Wilcoxon-Mann-Whitney comparison of frequency distributions")

  // Estimation
  .def("parametric_estimation", FrequencyDistributionWrap::parametric_estimation,
      return_value_policy<manage_new_object> (), "Parametric model estimation")
  .def("discrete_mixture_estimation1", FrequencyDistributionWrap::discrete_mixture_estimation1,
      return_value_policy<manage_new_object> (), "Discrete mixture estimation")
  .def("discrete_mixture_estimation2", FrequencyDistributionWrap::discrete_mixture_estimation2,
      return_value_policy<manage_new_object> (), "Discrete mixture estimation")
  .def("convolution_estimation1", FrequencyDistributionWrap::convolution_estimation1,
      return_value_policy<manage_new_object> (), "Convolution estimation")
  .def("convolution_estimation2", FrequencyDistributionWrap::convolution_estimation2,
      return_value_policy<manage_new_object> (), "Convolution estimation")
  .def("compound_estimation1", FrequencyDistributionWrap::compound_estimation1,
      return_value_policy<manage_new_object> (), "Compound distribution  estimation")
  .def("compound_estimation2", FrequencyDistributionWrap::compound_estimation2,
      return_value_policy<manage_new_object> (), "Compound distribution estimation")

  // Select
  .def("value_select", FrequencyDistributionWrap::value_select, return_value_policy<
      manage_new_object> (), args("min", "max", "keep"),
      "Selection of individuals according to the values taken by a variable")

  // Cluster
  .def("cluster_step", FrequencyDistributionWrap::cluster_step, return_value_policy<
      manage_new_object> (), args("step"), "See Cluster")
  .def("cluster_information", FrequencyDistributionWrap::cluster_information,
      return_value_policy<manage_new_object> (), args("ratio"),
      "Cluster with information")
  .def("cluster_limit", FrequencyDistributionWrap::cluster_limit, return_value_policy<
      manage_new_object> (), args("limits"), "See Cluster")
  .def("transcode", FrequencyDistributionWrap::transcode, return_value_policy<
      manage_new_object> (), "See Transcode")

  // Others
  .def("merge", FrequencyDistributionWrap::merge_frequency_distributions, return_value_policy<
      manage_new_object> (), args("frequency distributions"), "Merge frequency distributions")
  .def("shift", FrequencyDistributionWrap::shift,
      return_value_policy<manage_new_object> (), args("shift_param"),
      "Shift FrequencyDistribution")
  .def("fit", FrequencyDistributionWrap::fit, return_value_policy<manage_new_object> (),
      args("parametric_model"), "Fit frequency distribution")
  .def("get_plotable_list", FrequencyDistributionWrap::get_plotable, return_value_policy<
      manage_new_object> (), "Return a plotable for a list of frequency distributions")
  .def("get_plotable", &StatInterface::get_plotable, return_value_policy<
      manage_new_object> (), "Return a plotable (no parameters)")

  ;

  /*

   bool concentration_flag = false) const;
   bool plot_print(const char *path , int nb_histo = 0 ,  const FrequencyDistribution **histo = 0) const;
   bool plot_print(const char *path , double *cumul , double *concentration ,  double shift = 0.) const;
   bool survival_plot_print(const char *path , double *survivor) const;
   std::ostream& plot_title_print(std::ostream &os) const   { return os; }

   void plotable_frequency_write(SinglePlot &plot) const;
   void plotable_mass_write(SinglePlot &plot) const;
   void plotable_cumul_write(SinglePlot &plot , double *icumul = 0 ,   double scale = D_DEFAULT) const;
   void plotable_cumul_matching_write(SinglePlot &plot , int reference_offset ,   int  reference_nb_value , double *reference_cumul ,   double *icumul = 0) const;
   void plotable_concentration_write(SinglePlot &plot , double *icumul = 0 ,   double scale = D_DEFAULT) const;
   void plotable_survivor_write(SinglePlot &plot) const;

   double* cumul_computation(double scale = D_DEFAULT) const;
   double concentration_computation() const;

   double* survivor_function_computation(double scale = D_DEFAULT) const;
   double* concentration_function_computation(double scale = D_DEFAULT) const;

   Test* kruskal_wallis_test(int nb_histo , const FrequencyDistribution **ihisto) const;

   std::ostream& dissimilarity_ascii_write(std::ostream &os , int nb_histo ,  const FrequencyDistribution **ihisto ,   int type , double **dissimilarity) const;
   bool dissimilarity_ascii_write(StatError &error , const char *path ,   int nb_histo , const FrequencyDistribution **ihisto ,   int type , double **dissimilarity) const;
   bool dissimilarity_spreadsheet_write(StatError &error , const char *path ,   int nb_histo , const FrequencyDistribution **ihisto ,   int type , double **dissimilarity) const;

   void update(const Reestimation<double> *reestim , int inb_element);
   FrequencyDistribution* frequency_scale(int inb_element) const;

   FrequencyDistribution(int inb_value = 0)      :Reestimation<int>(inb_value) {}
   FrequencyDistribution(const Distribution &dist)      :Reestimation<int>(dist.nb_value) {}
   FrequencyDistribution(int inb_element , int *pelement);
   FrequencyDistribution(const FrequencyDistribution &histo , char transform , int param , int mode = FLOOR);


   bool comparison(StatError &error , std::ostream &os , int nb_histo ,
   const FrequencyDistribution **ihisto , int type , const char *path = 0 ,   char format = 'a') const;



   */

}

////////////////////// Export DiscreteDistributionData ////////////////////////////////

// Wrapper class

class DistributionDataWrap
{

public:

  // FrequencyDistribution constructor from a python list
  static boost::shared_ptr<DiscreteDistributionData>
  distribution_data_from_list(boost::python::list& int_list)
  {
    DiscreteDistributionData *histo = NULL;
    int nb_element = boost::python::len(int_list);

    if (!nb_element)
      stat_tool::wrap_util::throw_error("At least one observation"
        "is required to initialize FrequencyDistribution");

    stat_tool::wrap_util::auto_ptr_array<int> pelement(new int[nb_element]);
    for (int i = 0; i < nb_element; i++)
      pelement[i] = extract<int> (int_list[i]);

    // create new distribution_data object
    histo = new DiscreteDistributionData(nb_element, pelement.get());
    if (!histo)
      stat_tool::wrap_util::throw_error(
          "Could not initialize FrequencyDistribution from argument");

    return boost::shared_ptr<DiscreteDistributionData>(histo);
  }

  static boost::shared_ptr<DiscreteDistributionData>
  distribution_data_from_file(char* filename)
  {
    StatError error;
    DiscreteDistributionData *histo = NULL;
    histo = frequency_distribution_ascii_read(error, filename);

    if (!histo)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<DiscreteDistributionData>(histo);

  }

  WRAP_METHOD_FILE_ASCII_WRITE(DiscreteDistributionData);
  WRAP_METHOD_SURVIVAL_ASCII_WRITE(DiscreteDistributionData);
  WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(DiscreteDistributionData);
  WRAP_METHOD_SURVIVAL_PLOT_WRITE(DiscreteDistributionData);
  WRAP_METHOD0(DiscreteDistributionData, extract_model, DiscreteParametricModel);
  WRAP_METHOD0(DiscreteDistributionData, survival_get_plotable, MultiPlotSet);


};



// Boost declaration

void class_distribution_data()
{
    // _DiscreteDistributionData
  class_<DiscreteDistributionData, bases<FrequencyDistribution, StatInterface> >
    ("_DiscreteDistributionData", "_DiscreteDistributionData", init< boost::python::optional< int > >())

    .def("__init__", make_constructor(DistributionDataWrap::distribution_data_from_list ))
    .def("__init__", make_constructor(DistributionDataWrap::distribution_data_from_file))

    .def(self_ns::str(self)) // __str__

    // Output
    .def("survival_ascii_write", DistributionDataWrap::survival_ascii_write, "Return a string containing the object description (survival viewpoint)")
    .def("survival_plot_write", DistributionDataWrap::survival_plot_write, args("prefix", "title"), "Write GNUPLOT files (survival viewpoint)")
    .def("survival_spreadsheet_write", DistributionDataWrap::survival_spreadsheet_write, args("filename"), "Write object to filename (spreadsheet format)")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", DistributionDataWrap::survival_get_plotable, "Return a plotable object")
    .def("file_ascii_write", DistributionDataWrap::file_ascii_write, "Save frequency distribution into a file")

    DEF_RETURN_VALUE_NO_ARGS("extract_model", DistributionDataWrap::extract_model, "Return the 'model' part of the frequency distribution")

    ;

/*
    DiscreteDistributionData(const Distribution &dist)    :FrequencyDistribution(dist) { distribution = 0; }
    DiscreteDistributionData(const FrequencyDistribution &histo)    :FrequencyDistribution(histo) { distribution = 0; }
    DiscreteDistributionData(int inb_element , int *pelement)    :FrequencyDistribution(inb_element , pelement) { distribution = 0; }
    DiscreteDistributionData(const FrequencyDistribution &histo , char transform , int param , int mode = FLOOR)    :FrequencyDistribution(histo , transform , param , mode) { distribution = 0; }
    DiscreteDistributionData(int nb_histo , const FrequencyDistribution **phisto)    :FrequencyDistribution(nb_histo , phisto) { distribution = 0; }
    DiscreteDistributionData(const FrequencyDistribution &histo , const Distribution *dist);
    DiscreteDistributionData(const FrequencyDistribution &histo , const DiscreteParametric *dist);
    DiscreteDistributionData(const DiscreteDistributionData &histo , bool model_flag = true);
    DiscreteParametric* get_distribution() const { return distribution; }
*/

}
