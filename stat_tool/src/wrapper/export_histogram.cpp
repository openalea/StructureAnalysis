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
#include "export_base.h"

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


////////////////////// Export Histogram ////////////////////////////////////////

// Wrapper class class
class HistogramWrap
{

public:

  // __len__
  static int histo_get_nb_value(Histogram& h)
  {
    return h.nb_value;
  }

  static int histo_get_nb_element(Histogram& h)
  {
    return h.nb_element;
  }

  // __getitem__
  static int histo_get_item(Histogram& h, const object& key)
  {
    int val = extract<int>(key);

    if( (val >= 0) && (val < h.nb_value))
      return h.frequency[val];
    else
      {
	PyErr_SetString(PyExc_IndexError, "Bad Index");
	throw_error_already_set();
      }
  }


  static std::string f_comparison(const Histogram& h , const Histogram &histo)
  {
    ostringstream output;
    h.F_comparison(output, histo);
    return output.str();
      
  }

  static std::string t_comparison(const Histogram& h, const Histogram &histo)
  {
    ostringstream output;
    h.t_comparison(output, histo);
    return output.str();
  }

  static std::string wmw_comparison(const Histogram& h, const Histogram &histo)
  {
    ostringstream output;
    Format_error error;
    if(!h.wilcoxon_mann_whitney_comparison(error, output, histo))
      {
	stat_tool::wrap_util::throw_error(error);
      }
    return output.str();
      
  }

  static std::string comparison(const Histogram& h, 
				boost::python::tuple& histos, int type, 
				const char* filename, char format)
  {
    ostringstream output;
    Format_error error;
    int nb_histo = boost::python::len(histos);

    // Test list length
    if(nb_histo ==0)
	stat_tool::wrap_util::throw_error("No histogram to compare");

    stat_tool::wrap_util::auto_ptr_array<const Histogram*>
      ihisto(new const Histogram*[nb_histo]);

    // Extract list element
    for(int i=0; i<nb_histo; i++)
      ihisto[i] = boost::python::extract< Histogram *>(histos[i]);
    
    // call comparaison
    bool res = h.comparison(error, output, nb_histo, ihisto.get(), type, filename, format);
    

    if(!res)
      stat_tool::wrap_util::throw_error(error);

    return output.str();
  }


  // Estimation functions

  static Parametric_model* parametric_estimation(const Histogram& h, int ident ,
						 int min_inf_bound, bool flag)
  {
    Parametric_model* ret;
    Format_error error;
    
    ret = h.parametric_estimation(error, ident, min_inf_bound, flag);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Mixture* mixture_estimation_1(const Histogram& h, const Mixture& imixt, 
				       boost::python::list& estimate_tuple, 
				       int min_inf_bound,
				       bool flag, bool component_flag)
  {
    Mixture* ret;
    Format_error error;
    bool estimate[MIXTURE_NB_COMPONENT];

    //Parse estimate list
    int nb_element = boost::python::len(estimate_tuple);

    for (int i=0; i<nb_element; i++)
      estimate[i] = boost::python::extract< bool >(estimate_tuple[i]);
  
    ret = h.mixture_estimation(error, imixt, estimate, min_inf_bound, flag,
			     component_flag);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;

  }
  

  static Mixture* mixture_estimation_2(const Histogram& h, 
				       boost::python::list& ident_tuple, 
				       int min_inf_bound,
				       bool flag, bool component_flag, int penalty)
  {
    Mixture* ret;
    ostringstream output;
    Format_error error;

    int nb_component = boost::python::len(ident_tuple);
    int ident[MIXTURE_NB_COMPONENT]; 

    for (int i=0; i<nb_component; i++)
      ident[i] = boost::python::extract<int>(ident_tuple[i]);

    ret = h.mixture_estimation(error, output, 1, nb_component, ident, min_inf_bound,
			     flag, component_flag, penalty);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  
  static Convolution* convolution_estimation_1(const Histogram& h, 
					       const Parametric &known_dist ,
					       const Parametric &unknown_dist , 
					       int estimator, int nb_iter, 
					       double weight, int penalty_type, 
					       int outside)
  {
    Convolution* ret;
    ostringstream output;
    Format_error error;
    
    ret = h.convolution_estimation(error, output, known_dist, unknown_dist, estimator,
				 nb_iter, weight, penalty_type, outside);
    
    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Convolution* convolution_estimation_2(const Histogram& h, 
					       const Parametric &known_dist,
					       int min_inf_bound, 
					       int estimator, int nb_iter, 
					       double weight,  int penalty_type, 
					       int outside)
  {
    Convolution* ret;
    ostringstream output;
    Format_error error;
    
    ret = h.convolution_estimation(error, output, known_dist,  min_inf_bound, estimator,
				 nb_iter, weight, penalty_type, outside);
    
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }

  

  static Compound* compound_estimation_1(const Histogram& h, 
					 const Parametric &sum_dist,
					 const Parametric &dist,
					 char type,
					 int estimator, int nb_iter, 
					 double weight,  int penalty_type, 
					 int outside)
  {
    Compound* ret;
    ostringstream output;
    Format_error error;
    
    ret = h.compound_estimation(error, output, sum_dist,  dist, type, estimator,
				   nb_iter, weight, penalty_type, outside);
    
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }

  
  static Compound* compound_estimation_2(const Histogram& h, 
					 const Parametric &known_dist,
					 char type, int min_inf_bound,
					 int estimator, int nb_iter, 
					 double weight, int penalty_type, 
					 int outside)
  {
    Compound* ret;
    ostringstream output;
    Format_error error;
    
    ret = h.compound_estimation(error, output, known_dist, type, min_inf_bound, estimator,
			      nb_iter, weight, penalty_type, outside);
    
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }


  // Merge Histogram in a distribution_data
  static Distribution_data* merge_histograms(const Histogram& h, 
					     const boost::python::list& histo_list)
  {
    Distribution_data *histo = NULL;
    int nb_element = boost::python::len(histo_list) + 1;
    
    stat_tool::wrap_util::auto_ptr_array<const Histogram*>
      pelement(new const Histogram*[nb_element]);

    pelement[0] = &h;

    for (int i = 1; i < nb_element; i++)
      pelement[i] = extract<Histogram*>(histo_list[i-1]);
      
    
    // create new distribution_data object
    histo = new Distribution_data(nb_element, pelement.get());
    if (! histo)
      stat_tool::wrap_util::throw_error("Could not initialize Histogram from arguments");
	
    return histo;
  }


  static Distribution_data* value_select(const Histogram& h, 
					 int min, int max,
					 bool keep)
  {
    Format_error error;
    Distribution_data* ret = h.value_select(error, min, max, keep);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Distribution_data* shift(const Histogram& h, int shift_param)
  {
    Format_error error;
    Distribution_data* ret = h.shift(error, shift_param);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Parametric_model* fit(const Histogram& h, const Parametric& p)
  {
    Format_error error;
    Parametric_model* ret = h.fit(error, p);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  // Cluster
  static Distribution_data* cluster_step(const Histogram& h, int step)
  {
    Format_error error;
    Distribution_data* ret = h.cluster(error, step);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data* cluster_information(const Histogram& h, float ratio)
  {
    Format_error error;
    std::stringstream output;

    Distribution_data* ret = h.cluster(error, ratio, output);
    //print(output);

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }

  static Distribution_data* cluster_limit(const Histogram& h, 
					  boost::python::list& limit
					  )
  {

    Format_error error;

    int nb_limit = len(limit);
    stat_tool::wrap_util::auto_ptr_array<int> l (new int[nb_limit]);

    for (int i=0; i<nb_limit; i++)
	l[i] = extract<int>(limit[i]);
    
    Distribution_data* ret = h.cluster(error, nb_limit, l.get());


    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }


  static Distribution_data* transcode(const Histogram& h, 
				      boost::python::list& symbol
				      )
  {

    Format_error error;

    int nb_symbol = len(symbol);
    stat_tool::wrap_util::auto_ptr_array<int> l (new int[nb_symbol]);

    int expected_nb_symbol = h.nb_value - h.offset;
    if(nb_symbol != expected_nb_symbol)
      stat_tool::wrap_util::throw_error("Bad number of Symbol");

    for (int i=0; i<nb_symbol; i++)
	l[i] = extract<int>(symbol[i]);
    
    Distribution_data* ret = h.transcode(error, l.get());

    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    return ret;
  }



};



// Boost declaration

void class_histogram()
{

    enum_<stat_tool::wrap_util::UniqueInt<4, 1> >("VariableType")
      .value("ORDINAL", ORDINAL)
      .value("NUMERIC", NUMERIC)
      .value("SYMBOLIC", SYMBOLIC)
      .value("CIRCULAR", CIRCULAR)
      .export_values()
      ;
    
    enum_<stat_tool::wrap_util::UniqueInt<6, 2> >("LikelihoodPenaltyType")
      .value("AIC", AIC)
      .value("AICc", AICc)
      .value("BIC", BIC)
      .value("BICc", BICc)
      .value("ICL", ICL)
      .value("ICLc", ICLc)
      .export_values()
      ;

    enum_<stat_tool::wrap_util::UniqueInt<3, 3> >("SmoothingPenaltyType")
      .value("FIRST_DIFFERENCE", FIRST_DIFFERENCE)
      .value("SECOND_DIFFERENCE", SECOND_DIFFERENCE)
      .value("ENTROPY", ENTROPY)
      .export_values()
      ;
   
    enum_<stat_tool::wrap_util::UniqueInt<2, 4> >("OutsideType")
      .value("ZERO", ZERO)
      .value("CONTINUATION", CONTINUATION)
      .export_values()
      ;
   

    enum_<stat_tool::wrap_util::UniqueInt<4, 5> >("EstimatorType")
      .value("ZERO", ZERO)
      .value("LIKELIHOOD", LIKELIHOOD)
      .value("PENALIZED_LIKELIHOOD", PENALIZED_LIKELIHOOD)
      .value("PARAMETRIC_REGULARIZATION", PARAMETRIC_REGULARIZATION)
      .export_values()
      ;

    

    // _Histogram
    class_<Histogram>("_Histogram", no_init)
      .def( self == self )
      .def( self != self )
      .def(self_ns::str(self))
      .def("__len__", HistogramWrap::histo_get_nb_element)
      .def("__getitem__", HistogramWrap::histo_get_item)
      
      // Comparison
      .def("compare", HistogramWrap::comparison,
	   python::arg("histogram"),
	   "Comparison of histograms")

      .def("f_comparison", HistogramWrap::f_comparison,
	   python::arg("histogram"),
	   "F comparison of histograms")

      .def("t_comparison", HistogramWrap::t_comparison,
	   python::arg("histogram"),
	   "T comparison of histograms")

      .def("wmw_comparison", HistogramWrap::wmw_comparison,
	   python::arg("histogram"),
	   " Wilcoxon-Mann-Whitney comparison of histograms")
      
      // Estimation
      .def("parametric_estimation", HistogramWrap::parametric_estimation, 
	   return_value_policy< manage_new_object >(),
	   "Parametric model estimation")

      .def("mixture_estimation", HistogramWrap::mixture_estimation_1, 
	   return_value_policy< manage_new_object >(),
	   "Mixture Estimation")

      .def("mixture_estimation", HistogramWrap::mixture_estimation_2,
	   return_value_policy< manage_new_object >(),
	   "Mixture Estimation")

      .def("convolution_estimation", HistogramWrap::convolution_estimation_1, 
	   return_value_policy< manage_new_object >(),
	   "Convolution Estimation")

      .def("convolution_estimation", HistogramWrap::convolution_estimation_2,
	   return_value_policy< manage_new_object >(),
	   "Convolution Estimation")

      .def("compound_estimation", HistogramWrap::compound_estimation_1, 
	   return_value_policy< manage_new_object >(),
	   "Compound  Estimation")

      .def("compound_estimation", HistogramWrap::compound_estimation_2,
	   return_value_policy< manage_new_object >(),
	   "Compound  Estimation")

      // Select
      .def("value_select", HistogramWrap::value_select,
	   return_value_policy< manage_new_object >(),
	   python::args("min", "max", "keep"),
	   "Selection of individuals according to the values taken by a variable")

      // Cluster
      .def("cluster_step", HistogramWrap::cluster_step, 
	   return_value_policy< manage_new_object >(),
	   python::arg("step"),
	   "Cluster with step")

      .def("cluster_information", HistogramWrap::cluster_information, 
	   return_value_policy< manage_new_object >(),
	   python::arg("ratio"),
	   "Cluster with information")
      
      .def("cluster_limit", HistogramWrap::cluster_limit, 
	   return_value_policy< manage_new_object >(),
	   python::arg("limites"),
	   "Cluster with limits")

      .def("transcode", HistogramWrap::transcode,
	   return_value_policy< manage_new_object >(),
	   "Transcode")


      // Others
      .def("merge", HistogramWrap::merge_histograms,  
	   return_value_policy< manage_new_object >(),
	   python::arg("histograms"),
	   "Merge histograms")

      .def("shift", HistogramWrap::shift, 
	   return_value_policy< manage_new_object >(),
	   python::arg("shift_param"),
	   "Shift Histogram")

      .def("fit", HistogramWrap::fit, 
	   return_value_policy< manage_new_object >(),
	   python::arg("parametric_model"),
	   "Fit histogram")

      ;

    

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
    
    if (! nb_element)
      stat_tool::wrap_util::throw_error("At least one observation"
					"is required to initialize Histogram"); 
    
    stat_tool::wrap_util::auto_ptr_array<int> pelement(new int[nb_element]);
    for (int i = 0; i < nb_element; i++)
      pelement[i] = extract<int>(int_list[i]);
      
    
    // create new distribution_data object
    histo= new Distribution_data(nb_element, pelement.get());
    if (! histo)
      stat_tool::wrap_util::throw_error("Could not initialize Histogram from argument");
	
    return boost::shared_ptr<Distribution_data>(histo);
  }


  static boost::shared_ptr<Distribution_data> distribution_data_from_file(char* filename)
  {
    Format_error error;
    Distribution_data *histo = NULL;
    histo = histogram_ascii_read(error, filename);

    if(!histo)
	stat_tool::wrap_util::throw_error(error);


    return boost::shared_ptr<Distribution_data>(histo);

    }


  static std::string survival_ascii_write(const Distribution_data& d)
  {
    std::stringstream s;
    std::string res;
    
    d.survival_ascii_write(s);
    res = s.str();
    
    return res;
  }

  
  static void survival_spreadsheet_write(const Distribution_data& d,
						const std::string& filename)
  {
    Format_error error;

    if(!d.survival_spreadsheet_write(error, filename.c_str()))
      {
	stat_tool::wrap_util::throw_error(error);
      }

  }

   
  static void survival_plot_write(const Distribution_data& d,
				  const std::string& prefix, const std::string& title)
  {
    Format_error error;

    if(!d.survival_plot_write(error, prefix.c_str(), title.c_str()))
      {
	stat_tool::wrap_util::throw_error(error);
      }
  }


  static MultiPlotSet* survival_get_plotable(const Distribution_data& d)
  {
    Format_error error;
    MultiPlotSet* ret = d.survival_get_plotable(error);
    if(!ret)
      stat_tool::wrap_util::throw_error(error);
    
    return ret;
  }


  static Parametric_model* extract_model(const Distribution_data& d)
  {

    Format_error error;
    Parametric_model* ret = NULL;

    ret = d.extract_model(error);
    if(!ret) stat_tool::wrap_util::throw_error(error);

    return ret;
  }
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
    .def("survival_ascii_write", DistributionDataWrap::survival_ascii_write,
    	 "Return a string containing the object description (survival viewpoint)")

    .def("survival_plot_write", DistributionDataWrap::survival_plot_write,
	 python::args("prefix", "title"),
	 "Write GNUPLOT files (survival viewpoint)")

    .def("survival_get_plotable", DistributionDataWrap::survival_get_plotable,
	  return_value_policy< manage_new_object >(),
	 "Return a plotable object")

    .def("survival_spreadsheet_write", DistributionDataWrap::survival_spreadsheet_write,
	 python::arg("filename"),
	 "Write object to filename (spreadsheet format)")
    
    // Extract
    .def("extract_model", DistributionDataWrap::extract_model, 
	 return_value_policy< manage_new_object >(),
	 "Return the 'model' part of the histogram")

    ;


}
