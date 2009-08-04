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

#include "export_base.h"
#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"

// related to the enumerates
#include "stat_tool/regression.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distribution.h"
#include "stat_tool/distance_matrix.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"
using namespace boost::python;
using namespace boost;




// Overloads

// Format_error
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_update_overloads_1_3, update, 1, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_correction_update_overloads_2_4,
				       correction_update, 2, 4)


// Boost.Python Wrapper export function

void class_constant()
{
  //constant
  scope().attr("I_DEFAULT") = I_DEFAULT;
  scope().attr("D_DEFAULT") = D_DEFAULT;
  scope().attr("D_INF") = D_INF;
  scope().attr("MAX_DIFF_BOUND") = MAX_DIFF_BOUND;
  scope().attr("MAX_MEAN") = MAX_MEAN;

  // constants in stat_tool

  /*
  MAX_INF_BOUND
  MAX_DIFF_BOUND
  MAX_MEAN
  B_PROBABILITY
  B_THRESHOLD
  P_THRESHOLD
  NB_THRESHOLD
  SAMPLE_NB_VALUE_COEFF
  INF_BOUND_MARGIN
  SUP_BOUND_MARGIN
  POISSON_RATIO
  POISSON_RANGE
  NB_VALUE_COEFF
  MIN_RANGE
  MAX_SURFACE
  DIST_NB_ELEMENT
  CHI2_FREQUENCY
  MIN_T_VALUE
  MARGINAL_MAX_VALUE
  SKEWNESS_ROUNDNESS
  NB_ERROR
  LINE_NB_CHARACTER
  ASCII_NB_VALUE
  ASCII_SPACE
  ASCII_ROUNDNESS
  SPREADSHEET_ROUNDNESS
  DISPLAY_NB_INDIVIDUAL
  PLOT_NB_DISTRIBUTION
  PLOT_NB_HISTOGRAM
  PLOT_ROUNDNESS
  PLOT_SHIFT
  PLOT_MAX_SHIFT
  TIC_THRESHOLD
  PLOT_MASS_THRESHOLD
  YSCALE = 1.4;
  */

  // constants in regression
  //REGRESSION_NB_VECTOR
  //NEIGHBORHOOD

  // constants in mixture
  /*
    MIXTURE_NB_COMPONENT
    MIXTURE_PARAMETER
    MIN_WEIGHT_STEP
    MAX_WEIGHT_STEP
    MIXTURE_COEFF
    MIXTURE_LIKELIHOOD_DIFF
    MIXTURE_NB_ITER

   */
  //constant in reestimation
  scope().attr("CUMUL_THRESHOLD") = CUMUL_THRESHOLD;
  //BISECTION_RATIO_THRESHOLD
  //BISECTION_NB_ITER

  // constants in distance_matrix
  /*
  ASCII_NB_INDIVIDUAL
  PLOT_YMARGIN
  DISTANCE_ROUNDNESS
  GLOBAL_NB_ITER
  PARTITIONING_NB_ITER_1
  PARTITIONING_NB_ITER_2
  */

  // constants in curves
  /*
  MAX_FREQUENCY
  MAX_RANGE
  PLOT_NB_CURVE
  PLOT_MIN_FREQUENCY
    */

  //cosntants compound
  /*COMPOUND_THRESHOLD = 0.99999;  // seuil sur la fonction de repartition
  COMPOUND_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
  COMPOUND_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
  COMPOUND_NB_ITER = 10000;    // nombre maximum d'iterations EM
  COMPOUND_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
  COMPOUND_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
  COMPOUND_COEFF = 10;         // coefficient arrondi estimateur
*/
//constants convolution.h
/*
  CONVOLUTION_NB_DISTRIBUTION = 10;  // nombre maximum de lois elementaires
  CONVOLUTION_THRESHOLD = 0.9999;  // seuil sur la fonction de repartition
  CONVOLUTION_INIT_PROBABILITY = 0.001;  // seuil pour l'initialisation de la probabilite
  CONVOLUTION_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
  CONVOLUTION_NB_ITER = 10000;  // nombre maximum d'iterations EM
  CONVOLUTION_DIFFERENCE_WEIGHT = 0.5;  // poids par defaut de la penalisation
  CONVOLUTION_ENTROPY_WEIGHT = 0.1;  // poids par defaut de la penalisation (cas de l'entropie)
  CONVOLUTION_COEFF = 10;      // coefficient arrondi estimateur
*/
  // constants in vector
  scope().attr("VECTOR_NB_VARIABLE") = VECTOR_NB_VARIABLE;
  scope().attr("DISTANCE_NB_VECTOR") = DISTANCE_NB_VECTOR;
  scope().attr("CONTINGENCY_NB_VALUE") = CONTINGENCY_NB_VALUE;
  scope().attr("DISPLAY_CONTINGENCY_NB_VALUE") = DISPLAY_CONTINGENCY_NB_VALUE;
  scope().attr("VARIANCE_ANALYSIS_NB_VALUE") = VARIANCE_ANALYSIS_NB_VALUE;
  scope().attr("DISPLAY_CONDITIONAL_NB_VALUE") = DISPLAY_CONDITIONAL_NB_VALUE;
  scope().attr("PLOT_NB_VALUE") = PLOT_NB_VALUE;
  scope().attr("PLOT_RANGE_RATIO") = PLOT_RANGE_RATIO;
  scope().attr("NB_SYMBOL") = NB_SYMBOL;


  // distance_matrix
  enum_<stat_tool::wrap_util::UniqueInt<3, 1> > ("CriterionType")
  .value("NEAREST_NEIGHBOR", NEAREST_NEIGHBOR)
  .value("FARTHEST_NEIGHBOR",FARTHEST_NEIGHBOR)
  .value("AVERAGING", AVERAGING)
  .export_values();


  // from vectors.h
  enum_<stat_tool::wrap_util::UniqueInt<2, 2> >("DistanceType")
  .value("ABSOLUTE_VALUE", ABSOLUTE_VALUE)
  .value("QUADRATIC", QUADRATIC)
  .export_values()
  ;

  // from stat_tools.h
  enum_<stat_tool::wrap_util::UniqueInt<4, 3> > ("FittingType")
  .value("STANDARD_NORMAL", STANDARD_NORMAL)
  .value("CHI2" ,CHI2)
  .value("FISHER", FISHER)
  .value("STUDENT", STUDENT)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<5, 4> >("RoundType")
  .value("FLOOR", FLOOR)
  .value("ROUND", ROUND)
  .value("CEIL", CEIL)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<6, 5> >("DistributionIdentifierType")
  .value("NON_PARAMETRIC", NONPARAMETRIC)
  .value("BINOMIAL",BINOMIAL)
  .value("POISSON",POISSON)
  .value("NEGATIVE_BINOMIAL",NEGATIVE_BINOMIAL)
  .value("UNIFORM",UNIFORM)
  .value("MULTINOMIAL", MULTINOMIAL)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<4, 6> >("VariableType")
  .value("SYMBOLIC", SYMBOLIC)
  .value("ORDINAL", ORDINAL)
  .value("NUMERIC", NUMERIC)
  .value("CIRCULAR", CIRCULAR)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<3, 7> >("PearsonType")
  .value("PEARSON", PEARSON)
  .value("SPEARMAN", SPEARMAN)
  .value("KENDALL", KENDALL)
  .value("SPEARMAN2", SPEARMAN2)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<3, 8> >("SmoothingPenaltyType")
  .value("FIRST_DIFFERENCE", FIRST_DIFFERENCE)
  .value("SECOND_DIFFERENCE", SECOND_DIFFERENCE)
  .value("ENTROPY", ENTROPY)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<4, 9> >("EstimatorType")
  .value("ZERO", ZERO)
  .value("LIKELIHOOD", LIKELIHOOD)
  .value("PENALIZED_LIKELIHOOD", PENALIZED_LIKELIHOOD)
  .value("PARAMETRIC_REGULARIZATION", PARAMETRIC_REGULARIZATION)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<2, 10> >("OutsideType")
  .value("ZERO", ZERO)
  .value("CONTINUATION", CONTINUATION)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<3, 11> >("TOBEDONEType")
  .value("PARTIAL_LIKELIHOOD", PARTIAL_LIKELIHOOD)
  .value("COMPLETE_LIKELIHOOD", COMPLETE_LIKELIHOOD)
  .value("KAPLAN_MEIER", KAPLAN_MEIER)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<3, 12> >("ClusterType")
  .value("COMPUTED", COMPUTED)
  .value("ESTIMATED", ESTIMATED)
  .value("ONE_STEP_LATE", ONE_STEP_LATE)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<6, 13> >("LikelihoodPenaltyType")
  .value("AIC", AIC)
  .value("AICc", AICc)
  .value("BIC", BIC)
  .value("BICc", BICc)
  .value("ICL", ICL)
  .value("ICLc", ICLc)
  .export_values()
  ;

  enum_<stat_tool::wrap_util::UniqueInt<3, 13> > ("AlgoType")
  .value("AGGLOMERATIVE", AGGLOMERATIVE)
  .value("DIVISIVE", DIVISIVE)
  .value("ORDERING", ORDERING)
  .export_values();

  enum_<stat_tool::wrap_util::UniqueInt<6, 14> >("VariableTypeBis")
  .value("INT_VALUE", INT_VALUE)
  .value("REAL_VALUE", REAL_VALUE)
  .value("STATE", STATE)
  .value("OLD_INT_VALUE", OLD_INT_VALUE)
  .value("NB_INTERNODE", NB_INTERNODE)
  .value("AUXILIARY", AUXILIARY)
  .export_values()
  ;

  // regression
  enum_<stat_tool::wrap_util::UniqueInt<4, 15> >("RegressionType")
  .value("STAT_LINEAR", STAT_LINEAR)
  .value("STAT_LOGISTIC", STAT_LOGISTIC)
  .value("STAT_MONOMOLECULAR", STAT_MONOMOLECULAR)
  .value("STAT_NONPARAMETRIC", STAT_NONPARAMETRIC)
  .export_values()
  ;

}

// Format Error

void class_format_error()
{
  // _Format_error
  class_< Format_error >("_FormatError", init< optional< int > >())
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

}

class StatInterfaceWrap
{
public:

  WRAP_METHOD_ASCII_WRITE(STAT_interface);
  WRAP_METHOD_PLOT_WRITE(STAT_interface);
  WRAP_METHOD_SPREADSHEET_WRITE(STAT_interface);
};



void class_stat_interface()
{
  class_< STAT_interface, boost::noncopyable > ("_StatInterface", no_init)
    .def("ascii_write", &StatInterfaceWrap::ascii_write, args("exhaustive"),
    	"Return a string containing the object description (exhaustive or not)")

    .def("plot_write", &StatInterfaceWrap::plot_write,	args("prefix", "title"),
  		"Write GNUPLOT files (with prefix)")

    .def("spreadsheet_write", &StatInterfaceWrap::spreadsheet_write, args("filename"),
   		"Write object to filename (spreadsheet format)")

    .def("get_plotable", &STAT_interface::get_plotable, return_value_policy< manage_new_object >(),
  		"Return a plotable object" )
    ;

}


void class_forward()
{
  class_<Forward , boost::noncopyable ,bases<Parametric> > ("_Forward", no_init)

;
/*
    Forward(int inb_value = 0 , int iident = NONPARAMETRIC ,  int iinf_bound = I_DEFAULT , int isup_bound = I_DEFAULT , double iparameter = D_DEFAULT , double iprobability = D_DEFAULT)    :Parametric(inb_value , iident , iinf_bound , isup_bound , iparameter , iprobability) {}
    Forward(const Parametric &dist , int ialloc_nb_value = I_DEFAULT) :Parametric(dist , 'c' , ialloc_nb_value) { computation(dist); }
    Forward(const Forward &forward , int ialloc_nb_value = I_DEFAULT)  :Parametric((Parametric&)forward , 'c' , ialloc_nb_value) {}

    void computation(const Parametric &dist);
*/
}


