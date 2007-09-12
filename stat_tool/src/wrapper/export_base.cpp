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
 *        $Id: scenecontainer_wrap.cpp 559 2007-05-25 12:25:30Z dufourko $
 *                                                                       
 *-----------------------------------------------------------------------------*/

#include "export_base.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/convolution.h"
#include "stat_tool/compound.h"
#include "stat_tool/curves.h"
#include "stat_tool/mixture.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"


#include <boost/python.hpp>

using namespace boost::python;


// Wrappers



template <class U, class T, ostream & (T::*func)(ostream &, bool) const >
std::string ostream_converter(T &obj, bool exhaustive)
{  
   std::stringstream s;
   std::string res;
   Format_error error;

   obj.*func(s, exhaustive);
   res= s.str();

   return res;
}



// Overloads

// Format_error
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_update_overloads_1_3, update, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_correction_update_overloads_2_4, 
				       correction_update, 2, 4)


//histogram
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(value_select_overloads_3_4, value_select, 3, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(survival_plot_write_overloads_2_3, survival_plot_write, 2, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(comparison_overloads_5_7, comparison, 5, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(parametric_estimation_overloads_2_5, parametric_estimation, 2, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(type_parametric_estimation_overloads_1_4, type_parametric_estimation, 1, 4)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mixture_estimation_overloads_3_7, mixture_estimation, 3, 7)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mixture_estimation_overloads_2_6, mixture_estimation, 2, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(mixture_estimation_overloads_5_10, mixture_estimation, 5, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(convolution_estimation_overloads_4_9, convolution_estimation, 4, 9)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(compound_estimation_overloads_5_10, compound_estimation, 5, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(compound_estimation_overloads_4_10, compound_estimation, 4, 10)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(estimation_overloads_6_13, estimation, 6, 13)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(estimation_overloads_5_11, estimation, 5, 11)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(computation_overloads_0_2, computation, 0, 2)


// Boost.Python Wrapper export function
void class_base()
{
  // _Format_error
  class_< Format_error >("_Format_error", init< optional< int > >())
    .def("init", &Format_error::init)
    .def("update", &Format_error::update, Format_error_update_overloads_1_3())
    .def("correction_update", (void (Format_error::*)(const char*, const char*, int, int) )
	 &Format_error::correction_update, Format_error_correction_update_overloads_2_4())
    .def("correction_update", (void (Format_error::*)(const char*, int, int, int) )
	 &Format_error::correction_update, Format_error_correction_update_overloads_2_4())
    .def("get_nb_error", &Format_error::get_nb_error)
    .def("get_max_nb_error", &Format_error::get_max_nb_error)
    .def(self_ns::str(self))
    ;

  // Distribution base class
  //class_< Distribution, boost::noncopyable, Distribution_Wrapper >("Distribution", init< optional< int > >())
  class_< Distribution>("_Distribution", init< optional< int > >())
    .def(init< const Distribution&, double >())
    .def(init< const Histogram& >())
    .def(init< const Distribution&, optional< char, int > >())
    .def("survival_spreadsheet_write", &Distribution::survival_spreadsheet_write)
    .def("survival_plot_write", &Distribution::survival_plot_write, survival_plot_write_overloads_2_3())
    .def("mean_absolute_deviation_computation", &Distribution::mean_absolute_deviation_computation)
    .def("skewness_computation", &Distribution::skewness_computation)
    .def("kurtosis_computation", &Distribution::kurtosis_computation)
    .def("information_computation", &Distribution::information_computation)
    .def("first_difference_norm_computation", &Distribution::first_difference_norm_computation)
    .def("second_difference_norm_computation", &Distribution::second_difference_norm_computation)
    .def("likelihood_computation", (double (Distribution::*)(const Reestimation<int>&) const)&Distribution::likelihood_computation)
    .def("likelihood_computation", (double (Distribution::*)(const Reestimation<double>&) const)&Distribution::likelihood_computation)
    .def("chi2_fit", &Distribution::chi2_fit)
    .def("simulation", &Distribution::simulation)
    .def_readonly("offset", &Distribution::offset)
    .def_readonly("nb_value", &Distribution::nb_value)
    .def_readonly("max", &Distribution::max)
    .def_readonly("complement", &Distribution::complement)
    .def_readonly("mean", &Distribution::mean)
    .def_readonly("variance", &Distribution::variance)
    .def(self_ns::str(self))
    .def( self == self )
    .def( self != self )
    ;


  // Parametric base class
  class_< Parametric, bases< Distribution > >
    ("_Parametric", init< optional< int, int, int, int, double, double > >())
    .def(init< int, int, int, double, double, optional< double > >())
    .def(init< const Distribution&, optional< int > >())
    .def(init< const Distribution&, double >())
    //.def(init< const Parametric&, double >())
    .def(init< const Histogram& >())
    .def(init< const Parametric&, optional< char, int > >())
    .def("parametric_mean_computation", &Parametric::parametric_mean_computation)
    .def("parametric_variance_computation", &Parametric::parametric_variance_computation)
    .def("parametric_skewness_computation", &Parametric::parametric_skewness_computation)
    .def("parametric_kurtosis_computation", &Parametric::parametric_kurtosis_computation)
    .def("computation", &Parametric::computation, computation_overloads_0_2())
    .def("simulation", &Parametric::simulation)
    .def_readonly("ident", &Parametric::ident)
    .def_readonly("inf_bound", &Parametric::inf_bound)
    .def_readonly("sup_bound", &Parametric::sup_bound)
    .def_readonly("parameter", &Parametric::parameter)
    .def_readonly("probability", &Parametric::probability)
    .def(self_ns::str(self))
    ;


  
  // _Histogram
  class_<Histogram>("_Histogram", init<const Histogram&>())
    .def(init< optional< int > >())
    .def(init< const Distribution& >())
    .def(init< const Histogram&, char, int >())

    .def( self == self )
    .def( self != self )
    .def(self_ns::str(self))

  //   //Output
//     .def("ascii_write", 
// 	 (bool (Histogram::*)(Format_error&, const char*) const)&Histogram::ascii_write)
//     .def("ascii_write", &ostream_converter<Histogram, Histogram::ascii_write>)
//     .def("survival_spreadsheet_write", 
// 	 &Histogram::survival_spreadsheet_write)
//     .def("survival_spreadsheet_write", 
// 	 &Histogram::survival_spreadsheet_write)
//     .def("survival_plot_write", 
// 	 &Histogram::survival_plot_write, 
// 	 survival_plot_write_overloads_2_3())
    


    .def("shift", (Distribution_data* (Histogram::*)(Format_error&, int) const)&Histogram::shift, 
	 return_value_policy< manage_new_object >())

    .def("cluster", (Distribution_data* (Histogram::*)(Format_error&, int) const)&Histogram::cluster, 
	 return_value_policy< manage_new_object >())

    .def("cluster", 
	 (Distribution_data* (Histogram::*)(Format_error&, double, std::basic_ostream<char,std::char_traits<char> >&) const)&Histogram::cluster, 
	 return_value_policy< manage_new_object >())

    .def("cluster", (Distribution_data* (Histogram::*)(Format_error&, int, int*) const)&Histogram::cluster, 
	 return_value_policy< manage_new_object >())

    .def("transcode", &Histogram::transcode, return_value_policy< manage_new_object >())

    .def("value_select", &Histogram::value_select, 
	 return_value_policy< manage_new_object >(), value_select_overloads_3_4())

//     .def("build_time_events", &Histogram::build_time_events, 
// 	 return_value_policy< manage_new_object >())
    

    .def("comparison", &Histogram::comparison, comparison_overloads_5_7())

    .def("F_comparison", &Histogram::F_comparison)

    .def("t_comparison", &Histogram::t_comparison)

    .def("wilcoxon_mann_whitney_comparison", &Histogram::wilcoxon_mann_whitney_comparison)

    .def("fit", &Histogram::fit, return_value_policy< manage_new_object >())
    
    .def("parametric_estimation", 
	 (Parametric_model* (Histogram::*)(Format_error&, int, int, bool, double) const)&Histogram::parametric_estimation, 
	 return_value_policy< manage_new_object >(), 
	 parametric_estimation_overloads_2_5())
    
    .def("type_parametric_estimation", 
	 (Parametric_model* (Histogram::*)(Format_error&, int, bool, double) const)&Histogram::type_parametric_estimation, 
	 return_value_policy< manage_new_object >(), 
	 type_parametric_estimation_overloads_1_4())
    
    .def("mixture_estimation", 
	 (Mixture* (Histogram::*)(Format_error&, const Mixture&, bool*, int, bool, bool, double) const)&Histogram::mixture_estimation, 
	 return_value_policy< manage_new_object >(), 
	 mixture_estimation_overloads_3_7())
    
    .def("mixture_estimation", 
	 (Mixture* (Histogram::*)(Format_error&, const Mixture&, int, bool, bool, double) const)&Histogram::mixture_estimation, 
	 return_value_policy< manage_new_object >(), 
	 mixture_estimation_overloads_2_6())
    
    .def("mixture_estimation", 
	 (Mixture* (Histogram::*)(Format_error&, int, int*, int, bool, bool, double) const)&Histogram::mixture_estimation, 
	 return_value_policy< manage_new_object >(), 
	 mixture_estimation_overloads_3_7())
    
    .def("mixture_estimation", 
	 (Mixture* (Histogram::*)
	  (Format_error&, std::basic_ostream<char,std::char_traits<char> >&, int, int, int*, int, bool, bool, int, double)const
	  )&Histogram::mixture_estimation, 
	 return_value_policy< manage_new_object >(), 
	 mixture_estimation_overloads_5_10())

    .def("convolution_estimation", 
	 (Convolution* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, const Parametric&, int, int, double, int, int) const)&Histogram::convolution_estimation, 
	 return_value_policy< manage_new_object >(), 
	 convolution_estimation_overloads_4_9())

    .def("convolution_estimation", 
	 (Convolution* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, int, int, int, double, int, int) const)&Histogram::convolution_estimation, 
	 return_value_policy< manage_new_object >(), 
	 convolution_estimation_overloads_4_9())

    .def("compound_estimation", 
	 (Compound* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, const Parametric&, char, int, int, double, int, int) const)&Histogram::compound_estimation, 
	 return_value_policy< manage_new_object >(), 
	 compound_estimation_overloads_5_10())

    .def("compound_estimation", 
	 (Compound* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, char, int, int, int, double, int, int) const)&Histogram::compound_estimation, 
	 return_value_policy< manage_new_object >(), 
	 compound_estimation_overloads_4_10())

//     .def("estimation", 
// 	 (Parametric_model* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Histogram&, const Histogram&, const Histogram*, const Parametric&, int, int, int, double, int, int, double) const)&Histogram::estimation, 
// 	 return_value_policy< manage_new_object >(), estimation_overloads_6_13())

//     .def("estimation", 
// 	 (Parametric_model* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Histogram&, const Histogram&, const Histogram*, int, int, int, double, int, int) const)&Histogram::estimation, 
// 	 return_value_policy< manage_new_object >(), estimation_overloads_5_11())
    ;

}
