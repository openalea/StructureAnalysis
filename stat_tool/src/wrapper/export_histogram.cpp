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
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"
#include "export_base.h"

#include "stat_tool/stat_tools.h"

#include "stat_tool/distribution.h"
#include "stat_tool/convolution.h"
#include "stat_tool/mixture.h"
#include "stat_tool/compound.h"

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;

////////////////////// Export Histogram ////////////////////////////////////////


// Wrapper class class
class HistogramWrap
{

public:

  // __len__
  static int
  histo_get_nb_value(Histogram& h)
  {
    return h.nb_value;
  }

  static int
  histo_get_nb_element(Histogram& h)
  {
    return h.nb_element;
  }

  // __getitem__
  static int
  histo_get_item(Histogram& h, const object& key)
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
  f_comparison(const Histogram& h, const Histogram &histo)
  {
    ostringstream output;
    h.F_comparison(output, histo);
    return output.str();

  }

  static std::string
  t_comparison(const Histogram& h, const Histogram &histo)
  {
    ostringstream output;
    h.t_comparison(output, histo);
    return output.str();
  }

  static std::string
  wmw_comparison(const Histogram& h, const Histogram &histo)
  {
    ostringstream output;
    Format_error error;
    if (!h.wilcoxon_mann_whitney_comparison(error, output, histo))
      {
        stat_tool::wrap_util::throw_error(error);
      }
    return output.str();

  }

  static std::string
  comparison(const Histogram& h, boost::python::tuple& histos, int type,
      const char* filename, char format)
  {
    ostringstream output;
    Format_error error;
    int nb_histo = boost::python::len(histos);

    // Test list length
    if (nb_histo == 0)
      stat_tool::wrap_util::throw_error("No histogram to compare");

    stat_tool::wrap_util::auto_ptr_array<const Histogram*> ihisto(
        new const Histogram*[nb_histo]);

    // Extract list element
    for (int i = 0; i < nb_histo; i++)
      ihisto[i] = boost::python::extract<Histogram *>(histos[i]);

    // call comparaison
    bool res = h.comparison(error, output, nb_histo, ihisto.get(), type,
        filename, format);

    if (!res)
      stat_tool::wrap_util::throw_error(error);

    return output.str();
  }

  // Estimation functions

  static Mixture*
  mixture_estimation1(const Histogram& h, const Mixture& imixt,
      boost::python::list& estimate_tuple, int min_inf_bound, bool flag,
      bool component_flag)
  {
    Mixture* ret;
    Format_error error;
    bool estimate[MIXTURE_NB_COMPONENT];

    //Parse estimate list
    int nb_element = boost::python::len(estimate_tuple);

    for (int i = 0; i < nb_element; i++)
      estimate[i] = boost::python::extract<bool>(estimate_tuple[i]);

    ret = h.mixture_estimation(error, imixt, estimate, min_inf_bound, flag,
        component_flag);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;

  }

  static Mixture*
  mixture_estimation2(const Histogram& h, boost::python::list& ident_list,
      int min_inf_bound, bool flag, bool component_flag, int penalty)
  {
    Mixture* ret;
    ostringstream output;
    Format_error error;

    int nb_component = boost::python::len(ident_list);
    int ident[MIXTURE_NB_COMPONENT];

    for (int i = 0; i < nb_component; i++)
      ident[i] = boost::python::extract<int>(ident_list[i]);

    ret = h.mixture_estimation(error, output, 1, nb_component, ident,
        min_inf_bound, flag, component_flag, penalty);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    cout << output.str() << endl;

    return ret;
  }

  static Convolution*
  convolution_estimation1(const Histogram& h, const Parametric &known_dist,
      const Parametric &unknown_dist, int estimator, int nb_iter,
      double weight, int penalty_type, int outside)
  {
    Convolution* ret;
    ostringstream output;
    Format_error error;

    ret = h.convolution_estimation(error, output, known_dist, unknown_dist,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Convolution*
  convolution_estimation2(const Histogram& h, const Parametric &known_dist,
      int min_inf_bound, int estimator, int nb_iter, double weight,
      int penalty_type, int outside)
  {
    Convolution* ret;
    ostringstream output;
    Format_error error;

    ret = h.convolution_estimation(error, output, known_dist, min_inf_bound,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Compound*
  compound_estimation1(const Histogram& h, const Parametric &sum_dist,
      const Parametric &dist, char type, int estimator, int nb_iter,
      double weight, int penalty_type, int outside)
  {
    Compound* ret;
    ostringstream output;
    Format_error error;

    ret = h.compound_estimation(error, output, sum_dist, dist, type, estimator,
        nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Compound*
  compound_estimation2(const Histogram& h, const Parametric &known_dist,
      char type, int min_inf_bound, int estimator, int nb_iter, double weight,
      int penalty_type, int outside)
  {
    Compound* ret;
    ostringstream output;
    Format_error error;

    ret = h.compound_estimation(error, output, known_dist, type, min_inf_bound,
        estimator, nb_iter, weight, penalty_type, outside);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Merge Histogram in a distribution_data
  static Distribution_data*
  merge_histograms(const Histogram& h, const boost::python::list& histo_list)
  {
    Distribution_data *histo = NULL;
    int nb_element = boost::python::len(histo_list) + 1;

    stat_tool::wrap_util::auto_ptr_array<const Histogram*> pelement(
        new const Histogram*[nb_element]);

    pelement[0] = &h;

    for (int i = 1; i < nb_element; i++)
      pelement[i] = extract<Histogram*> (histo_list[i - 1]);

    // create new distribution_data object
    histo = new Distribution_data(nb_element, pelement.get());
    if (!histo)
      stat_tool::wrap_util::throw_error(
          "Could not initialize Histogram from arguments");

    return histo;
  }

  WRAP_METHOD3(Histogram, parametric_estimation, Parametric_model, int, int, bool);
  WRAP_METHOD3(Histogram, value_select, Distribution_data, int, int, bool);
  WRAP_METHOD1(Histogram, shift, Distribution_data, int);
  WRAP_METHOD1(Histogram, fit, Parametric_model, Parametric);

  // Cluster
  static Distribution_data*
  cluster_step(const Histogram& h, int step)
  {
    Format_error error;
    Distribution_data* ret = h.cluster(error, step);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data*
  cluster_information(const Histogram& h, float ratio)
  {
    Format_error error;
    std::stringstream output;

    Distribution_data* ret = h.cluster(error, ratio, output);

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data*
  cluster_limit(const Histogram& h, boost::python::list& limit)
  {

    Format_error error;

    int nb_limit = len(limit);
    stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_limit]);

    for (int i = 0; i < nb_limit; i++)
      l[i] = extract<int> (limit[i]);

    Distribution_data* ret = h.cluster(error, nb_limit + 1, l.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data*
  transcode(const Histogram& h, boost::python::list& symbol)
  {

    Format_error error;

    int nb_symbol = len(symbol);
    stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

    int expected_nb_symbol = h.nb_value - h.offset;
    if (nb_symbol != expected_nb_symbol)
      stat_tool::wrap_util::throw_error("Bad number of Symbol");

    for (int i = 0; i < nb_symbol; i++)
      l[i] = extract<int> (symbol[i]);

    Distribution_data* ret = h.transcode(error, l.get());

    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static MultiPlotSet*
  get_plotable(const Histogram& p, const boost::python::list& hist_list)
  {
    Format_error error;
    int nb_hist = boost::python::len(hist_list);
    stat_tool::wrap_util::auto_ptr_array<const Histogram *> hists(
        new const Histogram*[nb_hist]);

    for (int i = 0; i < nb_hist; i++)
      hists[i] = extract<const Histogram*> (hist_list[i]);

    const Histogram** d = hists.get();

    MultiPlotSet* ret = p.get_plotable_histograms(error, nb_hist, d);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

};



// Boost declaration

void class_histogram()
{

   // _Histogram
   class_<Histogram> ("_Histogram", no_init)
  .def(self == self)
  .def(self != self) .def(self_ns::str(self))
  .def("__len__", HistogramWrap::histo_get_nb_element)
  .def("__getitem__", HistogramWrap::histo_get_item)

  // Comparison
  .def("compare", HistogramWrap::comparison,
      args("histogram"),"Comparison of histograms")
  .def("f_comparison", HistogramWrap::f_comparison,
      args("histogram"), "F comparison of histograms")
  .def("t_comparison", HistogramWrap::t_comparison, args("histogram"),
      "T comparison of histograms")
  .def("wmw_comparison", HistogramWrap::wmw_comparison,
      args("histogram"),
      " Wilcoxon-Mann-Whitney comparison of histograms")

  // Estimation
  .def("parametric_estimation", HistogramWrap::parametric_estimation,
      return_value_policy<manage_new_object> (), "Parametric model estimation")
  .def("mixture_estimation1", HistogramWrap::mixture_estimation1,
      return_value_policy<manage_new_object> (), "Mixture Estimation")
  .def("mixture_estimation2", HistogramWrap::mixture_estimation2,
      return_value_policy<manage_new_object> (), "Mixture Estimation")
  .def("convolution_estimation1", HistogramWrap::convolution_estimation1,
      return_value_policy<manage_new_object> (), "Convolution Estimation")
  .def("convolution_estimation2", HistogramWrap::convolution_estimation2,
      return_value_policy<manage_new_object> (), "Convolution Estimation")
  .def("compound_estimation1", HistogramWrap::compound_estimation1,
      return_value_policy<manage_new_object> (), "Compound  Estimation")
  .def("compound_estimation2", HistogramWrap::compound_estimation2,
      return_value_policy<manage_new_object> (), "Compound  Estimation")

  // Select
  .def("value_select", HistogramWrap::value_select, return_value_policy<
      manage_new_object> (), args("min", "max", "keep"),
      "Selection of individuals according to the values taken by a variable")

  // Cluster
  .def("cluster_step", HistogramWrap::cluster_step, return_value_policy<
      manage_new_object> (), args("step"), "See Cluster")
  .def("cluster_information", HistogramWrap::cluster_information,
      return_value_policy<manage_new_object> (), args("ratio"),
      "Cluster with information")
  .def("cluster_limit", HistogramWrap::cluster_limit, return_value_policy<
      manage_new_object> (), args("limits"), "See Cluster")
  .def("transcode", HistogramWrap::transcode, return_value_policy<
      manage_new_object> (), "See Transcode")

  // Others
  .def("merge", HistogramWrap::merge_histograms, return_value_policy<
      manage_new_object> (), args("histograms"), "Merge histograms")
  .def("shift", HistogramWrap::shift,
      return_value_policy<manage_new_object> (), args("shift_param"),
      "Shift Histogram")
  .def("fit", HistogramWrap::fit, return_value_policy<manage_new_object> (),
      args("parametric_model"), "Fit histogram")
  .def("get_plotable_list", HistogramWrap::get_plotable, return_value_policy<
      manage_new_object> (), "Return a plotable for a list of histograms")
  .def("get_plotable", &STAT_interface::get_plotable, return_value_policy<
      manage_new_object> (), "Return a plotable (no parameters)")

  ;

  /*

   bool concentration_flag = false) const;
   bool plot_print(const char *path , int nb_histo = 0 ,  const Histogram **histo = 0) const;
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

   Test* kruskal_wallis_test(int nb_histo , const Histogram **ihisto) const;

   std::ostream& dissimilarity_ascii_write(std::ostream &os , int nb_histo ,  const Histogram **ihisto ,   int type , double **dissimilarity) const;
   bool dissimilarity_ascii_write(Format_error &error , const char *path ,   int nb_histo , const Histogram **ihisto ,   int type , double **dissimilarity) const;
   bool dissimilarity_spreadsheet_write(Format_error &error , const char *path ,   int nb_histo , const Histogram **ihisto ,   int type , double **dissimilarity) const;

   void update(const Reestimation<double> *reestim , int inb_element);
   Histogram* frequency_scale(int inb_element) const;

   Histogram(int inb_value = 0)      :Reestimation<int>(inb_value) {}
   Histogram(const Distribution &dist)      :Reestimation<int>(dist.nb_value) {}
   Histogram(int inb_element , int *pelement);
   Histogram(const Histogram &histo , char transform , int param , int mode = FLOOR);


   bool comparison(Format_error &error , std::ostream &os , int nb_histo ,
   const Histogram **ihisto , int type , const char *path = 0 ,   char format = 'a') const;



   */

}

////////////////////// Export Distribution_data ////////////////////////////////

// Wrapper class

class DistributionDataWrap
{

public:

  // Histogram constructor from a python list
  static boost::shared_ptr<Distribution_data>
  distribution_data_from_list(boost::python::list& int_list)
  {
    Distribution_data *histo = NULL;
    int nb_element = boost::python::len(int_list);

    if (!nb_element)
      stat_tool::wrap_util::throw_error("At least one observation"
        "is required to initialize Histogram");

    stat_tool::wrap_util::auto_ptr_array<int> pelement(new int[nb_element]);
    for (int i = 0; i < nb_element; i++)
      pelement[i] = extract<int> (int_list[i]);

    // create new distribution_data object
    histo = new Distribution_data(nb_element, pelement.get());
    if (!histo)
      stat_tool::wrap_util::throw_error(
          "Could not initialize Histogram from argument");

    return boost::shared_ptr<Distribution_data>(histo);
  }

  static boost::shared_ptr<Distribution_data>
  distribution_data_from_file(char* filename)
  {
    Format_error error;
    Distribution_data *histo = NULL;
    histo = histogram_ascii_read(error, filename);

    if (!histo)
      stat_tool::wrap_util::throw_error(error);

    return boost::shared_ptr<Distribution_data>(histo);

  }

  WRAP_METHOD_FILE_ASCII_WRITE(Distribution_data);
  WRAP_METHOD_SURVIVAL_ASCII_WRITE(Distribution_data);
  WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(Distribution_data);
  WRAP_METHOD_SURVIVAL_PLOT_WRITE(Distribution_data);
  WRAP_METHOD0(Distribution_data, extract_model, Parametric_model);
  WRAP_METHOD0(Distribution_data, survival_get_plotable, MultiPlotSet);


};



// Boost declaration

void class_distribution_data()
{
    // _Distribution_data
  class_<Distribution_data, bases<Histogram, STAT_interface> >
    ("_DistributionData", "_DistributionData", init< optional< int > >())

    .def("__init__", make_constructor(DistributionDataWrap::distribution_data_from_list ))
    .def("__init__", make_constructor(DistributionDataWrap::distribution_data_from_file))

    .def(self_ns::str(self)) // __str__

    // Output
    .def("survival_ascii_write", DistributionDataWrap::survival_ascii_write, "Return a string containing the object description (survival viewpoint)")
    .def("survival_plot_write", DistributionDataWrap::survival_plot_write, args("prefix", "title"), "Write GNUPLOT files (survival viewpoint)")
    .def("survival_spreadsheet_write", DistributionDataWrap::survival_spreadsheet_write, args("filename"), "Write object to filename (spreadsheet format)")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", DistributionDataWrap::survival_get_plotable, "Return a plotable object")
    .def("file_ascii_write", DistributionDataWrap::file_ascii_write, "Save histogram into a file")

    DEF_RETURN_VALUE_NO_ARGS("extract_model", DistributionDataWrap::extract_model, "Return the 'model' part of the histogram")

    ;

/*
    Distribution_data(const Distribution &dist)    :Histogram(dist) { distribution = 0; }
    Distribution_data(const Histogram &histo)    :Histogram(histo) { distribution = 0; }
    Distribution_data(int inb_element , int *pelement)    :Histogram(inb_element , pelement) { distribution = 0; }
    Distribution_data(const Histogram &histo , char transform , int param , int mode = FLOOR)    :Histogram(histo , transform , param , mode) { distribution = 0; }
    Distribution_data(int nb_histo , const Histogram **phisto)    :Histogram(nb_histo , phisto) { distribution = 0; }
    Distribution_data(const Histogram &histo , const Distribution *dist);
    Distribution_data(const Histogram &histo , const Parametric *dist);
    Distribution_data(const Distribution_data &histo , bool model_flag = true);
    Parametric* get_distribution() const { return distribution; }
*/

}
