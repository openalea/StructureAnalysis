

// Includes ====================================================================
#include "stat_tool/stat_tools.h"
#include "stat_tool/convolution.h"
#include "stat_tool/compound.h"
#include "stat_tool/curves.h"
#include "stat_tool/mixture.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distribution.h"

#include <boost/python.hpp>
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::len
#include <boost/python/make_constructor.hpp>
// definition of boost::python::make_constructor

// Using =======================================================================
using namespace boost::python;

// Declarations ================================================================
namespace  {

template<int num, int id> struct UniqueInt
{
   int v;
   enum { value=num };

   UniqueInt(int _v) : v(_v) { }
   operator int() const { return v; }
};

int ERROR_LENGTH_() { return ERROR_LENGTH; }
int I_DEFAULT_() { return I_DEFAULT; }
double D_DEFAULT_() { return D_DEFAULT; }
double D_INF_() { return D_INF; }
double DOUBLE_ERROR_() { return DOUBLE_ERROR; }
double SELF_TRANSITION_() { return SELF_TRANSITION; }

struct Distribution_Wrapper: Distribution
{
    Distribution_Wrapper(PyObject* self_):
        Distribution(), self(self_) {}

    Distribution_Wrapper(PyObject* self_, int p0):
        Distribution(p0), self(self_) {}

    Distribution_Wrapper(PyObject* self_, const Distribution& p0, double p1):
        Distribution(p0, p1), self(self_) {}

    Distribution_Wrapper(PyObject* self_, const Histogram& p0):
        Distribution(p0), self(self_) {}

    Distribution_Wrapper(PyObject* self_, const Distribution& p0):
        Distribution(p0), self(self_) {}

    Distribution_Wrapper(PyObject* self_, const Distribution& p0, char p1):
        Distribution(p0, p1), self(self_) {}

    Distribution_Wrapper(PyObject* self_, const Distribution& p0, char p1, int p2):
        Distribution(p0, p1, p2), self(self_) {}

    std::basic_ostream<char,std::char_traits<char> >& plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return call_method< std::basic_ostream<char,std::char_traits<char> >& >(self, "plot_title_print", p0);
    }

    std::basic_ostream<char,std::char_traits<char> >& default_plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return Distribution::plot_title_print(p0);
    }

    PyObject* self;
};

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Distribution_survival_plot_write_overloads_2_3, survival_plot_write, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_update_overloads_1_3, update, 1, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Format_error_correction_update_overloads_2_4, correction_update, 2, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_value_select_overloads_3_4, value_select, 3, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_survival_plot_write_overloads_2_3, survival_plot_write, 2, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_comparison_overloads_5_7, comparison, 5, 7)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_parametric_estimation_overloads_2_5, parametric_estimation, 2, 5)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_type_parametric_estimation_overloads_1_4, type_parametric_estimation, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_mixture_estimation_overloads_3_7, mixture_estimation, 3, 7)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_mixture_estimation_overloads_2_6, mixture_estimation, 2, 6)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_mixture_estimation_overloads_5_10, mixture_estimation, 5, 10)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_convolution_estimation_overloads_4_9, convolution_estimation, 4, 9)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_compound_estimation_overloads_5_10, compound_estimation, 5, 10)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_compound_estimation_overloads_4_10, compound_estimation, 4, 10)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_estimation_overloads_6_13, estimation, 6, 13)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Histogram_estimation_overloads_5_11, estimation, 5, 11)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Reestimation_int_parametric_estimation_overloads_1_4, parametric_estimation, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Reestimation_int_type_parametric_estimation_overloads_1_4, type_parametric_estimation, 1, 4)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Reestimation_int_type_parametric_estimation_overloads_0_3, type_parametric_estimation, 0, 3)

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Reestimation_int_state_occupancy_estimation_overloads_4_5, state_occupancy_estimation, 4, 5)

struct Parametric_Wrapper: Parametric
{
    Parametric_Wrapper(PyObject* self_):
        Parametric(), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0):
        Parametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1):
        Parametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2):
        Parametric(p0, p1, p2), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3):
        Parametric(p0, p1, p2, p3), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3, double p4):
        Parametric(p0, p1, p2, p3, p4), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, int p3, double p4, double p5):
        Parametric(p0, p1, p2, p3, p4, p5), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, double p3, double p4):
        Parametric(p0, p1, p2, p3, p4), self(self_) {}

    Parametric_Wrapper(PyObject* self_, int p0, int p1, int p2, double p3, double p4, double p5):
        Parametric(p0, p1, p2, p3, p4, p5), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0):
        Parametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0, int p1):
        Parametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Distribution& p0, double p1):
        Parametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Parametric& p0, double p1):
        Parametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Forward& p0):
        Parametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Parametric& p0):
        Parametric(p0), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Parametric& p0, char p1):
        Parametric(p0, p1), self(self_) {}

    Parametric_Wrapper(PyObject* self_, const Parametric& p0, char p1, int p2):
        Parametric(p0, p1, p2), self(self_) {}

    std::basic_ostream<char,std::char_traits<char> >& plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return call_method< std::basic_ostream<char,std::char_traits<char> >& >(self, "plot_title_print", p0);
    }

    std::basic_ostream<char,std::char_traits<char> >& default_plot_title_print(std::basic_ostream<char,std::char_traits<char> >& p0) const {
        return Parametric::plot_title_print(p0);
    }

    PyObject* self;
};


BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Parametric_computation_overloads_0_2, computation, 0, 2)

/*************************************************************
 *
 *  Wrappers:
 */


Distribution_data* Histogram_wrapper_init(boost::python::list int_list)
{
   int nb_element, i;
   int *pelement= NULL;
   object o;
   ostringstream error_message;
   bool status= true;
   Distribution_data *histo= NULL;
   // string s;

   nb_element= boost::python::len(int_list);
   if (nb_element > 0)
   {
      pelement= new int[nb_element];
      for (i= 0; i < nb_element; i++)
      {
         o= int_list[i];
         try
         {
            extract<int> x(o);
            if (x.check())
               pelement[i]= x();
            else
               status=false;
         }
         catch (...)
         {
            status= false;
         }
         if (!status)
            error_message << "incorrect type for element " << i
                          << " of argument list: expecting an integer" << endl;
      }
      if (!status)
      {
         delete [] pelement;
         pelement= NULL;
         PyErr_SetString(PyExc_TypeError, (error_message.str()).c_str());
         throw_error_already_set();
      }
      else
      {
         histo= new Distribution_data(nb_element, pelement);
         delete [] pelement;
         pelement= NULL;
         if (histo == NULL)
         {
            status= false;
            error_message << "could not initialize a Histogram object from argument"; // << endl;
            PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
            throw_error_already_set();
         }
      }
   }
   else
   {
      status= false;
      error_message << "at least one observation required to initialize a Histogram object"; // << endl;
      PyErr_SetString(PyExc_IndexError, (error_message.str()).c_str());
      throw_error_already_set();
   }
   return histo;
}

std::string Histogram_wrapper_display(const Distribution_data& histo,
                                      bool exhaustive)
{
   std::stringstream s;
   std::string res;
   Format_error error;

   histo.ascii_write(s, exhaustive);
   res= s.str();

   return res;
}

std::string Histogram_wrapper_display_survival(const Distribution_data& histo)
{
   std::stringstream s;
   std::string res;
   Format_error error;

   histo.survival_ascii_write(s);
   res= s.str();

   return res;
}

void Histogram_wrapper_ascii_write_file(const Distribution_data& histo,
                                        const char* path,
                                        bool exhaustive)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;
   ofstream out_file(path);

   // status= histo.ascii_write(error, path, exhaustive);
   // obviously the commented solution would be more satisfying,
   // but the "file" option is overridden for Parametric_model::ascii_write
   status= histo.Histogram::ascii_write(out_file , exhaustive , true);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

void Histogram_wrapper_spreadsheet_write(const Distribution_data& histo,
                                         const char *path)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= histo.spreadsheet_write(error, path);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

void Histogram_wrapper_plot_write(const Distribution_data& histo,
                                  const char* prefix,
                                  const char* title)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= histo.plot_write(error, prefix, title);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

void Histogram_wrapper_plot_write_survival(const Distribution_data& histo,
                                           const char* prefix,
                                           const char* title)
{
   bool status= true;
   ostringstream error_message;
   Format_error error;

   status= histo.survival_plot_write(error, prefix, title);
   if (not status)
   {
      error_message << error;
      PyErr_SetString(PyExc_RuntimeError, (error_message.str()).c_str());
      throw_error_already_set();
   }
}

}// namespace


// Module ======================================================================
BOOST_PYTHON_MODULE(cstat_tool)
{
    class_< Distribution, boost::noncopyable, Distribution_Wrapper >("Distribution", init< optional< int > >())
    // class_< Distribution>("Distribution", init< optional< int > >())
        .def(init< const Distribution&, double >())
        .def(init< const Histogram& >())
        .def(init< const Distribution&, optional< char, int > >())
        .def("survival_spreadsheet_write", &Distribution::survival_spreadsheet_write)
        .def("survival_plot_write", &Distribution::survival_plot_write, Distribution_survival_plot_write_overloads_2_3())
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

    class_< Format_error >("Format_error", init< optional< int > >())
        .def("init", &Format_error::init)
        .def("update", &Format_error::update, Format_error_update_overloads_1_3())
        .def("correction_update", (void (Format_error::*)(const char*, const char*, int, int) )&Format_error::correction_update, Format_error_correction_update_overloads_2_4())
        .def("correction_update", (void (Format_error::*)(const char*, int, int, int) )&Format_error::correction_update, Format_error_correction_update_overloads_2_4())
        .def("get_nb_error", &Format_error::get_nb_error)
        .def("get_max_nb_error", &Format_error::get_max_nb_error)
        .def(self_ns::str(self))
    ;

    class_< Histogram >("CHistogram", init< const Histogram& >());

    // class_< Distribution_data >("Histogram", init< const Histogram& >())
    class_< Distribution_data, bases< Histogram> >("Histogram", init< const Histogram& >())
        .def(init< optional< int > >())
        .def(init< const Distribution& >())
        .def(init< const Distribution_data& >())
        .def(init< int, int* >())
        .def(init< int, const Histogram** >())
        .def(init< const Histogram&, char, int >())
        .def("__init__", make_constructor(Histogram_wrapper_init))
        .def("shift", (Distribution_data* (Histogram::*)(Format_error&, int) const)&Histogram::shift, return_value_policy< manage_new_object >())
        .def("cluster", (Distribution_data* (Histogram::*)(Format_error&, int) const)&Histogram::cluster, return_value_policy< manage_new_object >())
        .def("cluster", (Distribution_data* (Histogram::*)(Format_error&, double, std::basic_ostream<char,std::char_traits<char> >&) const)&Histogram::cluster, return_value_policy< manage_new_object >())
        .def("cluster", (Distribution_data* (Histogram::*)(Format_error&, int, int*) const)&Histogram::cluster, return_value_policy< manage_new_object >())
        .def("transcode", &Histogram::transcode, return_value_policy< manage_new_object >())
        .def("value_select", &Histogram::value_select, return_value_policy< manage_new_object >(), Histogram_value_select_overloads_3_4())
        // .def("build_time_events", &Histogram::build_time_events, return_value_policy< manage_new_object >())
        .def("ascii_write", (bool (Histogram::*)(Format_error&, const char*) const)&Histogram::ascii_write)
        .def("survival_spreadsheet_write", &Histogram::survival_spreadsheet_write)
        .def("survival_plot_write", &Histogram::survival_plot_write, Histogram_survival_plot_write_overloads_2_3())
        .def("comparison", &Histogram::comparison, Histogram_comparison_overloads_5_7())
        .def("F_comparison", &Histogram::F_comparison)
        .def("t_comparison", &Histogram::t_comparison)
        .def("wilcoxon_mann_whitney_comparison", &Histogram::wilcoxon_mann_whitney_comparison)
        .def("fit", &Histogram::fit, return_value_policy< manage_new_object >())
        // .def("parametric_estimation", (Parametric_model* (Histogram::*)(Format_error&, int, int, bool, double) const)&Histogram::parametric_estimation, return_value_policy< manage_new_object >(), Histogram_parametric_estimation_overloads_2_5())
        // .def("type_parametric_estimation", (Parametric_model* (Histogram::*)(Format_error&, int, bool, double) const)&Histogram::type_parametric_estimation, return_value_policy< manage_new_object >(), Histogram_type_parametric_estimation_overloads_1_4())
        // .def("mixture_estimation", (Mixture* (Histogram::*)(Format_error&, const Mixture&, bool*, int, bool, bool, double) const)&Histogram::mixture_estimation, return_value_policy< manage_new_object >(), Histogram_mixture_estimation_overloads_3_7())
        // .def("mixture_estimation", (Mixture* (Histogram::*)(Format_error&, const Mixture&, int, bool, bool, double) const)&Histogram::mixture_estimation, return_value_policy< manage_new_object >(), Histogram_mixture_estimation_overloads_2_6())
        // .def("mixture_estimation", (Mixture* (Histogram::*)(Format_error&, int, int*, int, bool, bool, double) const)&Histogram::mixture_estimation, return_value_policy< manage_new_object >(), Histogram_mixture_estimation_overloads_3_7())
        // .def("mixture_estimation", (Mixture* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, int, int, int*, int, bool, bool, int, double) const)&Histogram::mixture_estimation, return_value_policy< manage_new_object >(), Histogram_mixture_estimation_overloads_5_10())
        // .def("convolution_estimation", (Convolution* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, const Parametric&, int, int, double, int, int) const)&Histogram::convolution_estimation, return_value_policy< manage_new_object >(), Histogram_convolution_estimation_overloads_4_9())
        // .def("convolution_estimation", (Convolution* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, int, int, int, double, int, int) const)&Histogram::convolution_estimation, return_value_policy< manage_new_object >(), Histogram_convolution_estimation_overloads_4_9())
        // .def("compound_estimation", (Compound* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, const Parametric&, char, int, int, double, int, int) const)&Histogram::compound_estimation, return_value_policy< manage_new_object >(), Histogram_compound_estimation_overloads_5_10())
        // .def("compound_estimation", (Compound* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Parametric&, char, int, int, int, double, int, int) const)&Histogram::compound_estimation, return_value_policy< manage_new_object >(), Histogram_compound_estimation_overloads_4_10())
        // .def("estimation", (Parametric_model* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Histogram&, const Histogram&, const Histogram*, const Parametric&, int, int, int, double, int, int, double) const)&Histogram::estimation, return_value_policy< manage_new_object >(), Histogram_estimation_overloads_6_13())
        // .def("estimation", (Parametric_model* (Histogram::*)(Format_error&, std::basic_ostream<char,std::char_traits<char> >&, const Histogram&, const Histogram&, const Histogram*, int, int, int, double, int, int) const)&Histogram::estimation, return_value_policy< manage_new_object >(), Histogram_estimation_overloads_5_11())
        // .def("state_occupancy_estimation", &Reestimation<int>::state_occupancy_estimation, Reestimation_int_state_occupancy_estimation_overloads_4_5())
        .def("ascii_write", &Histogram_wrapper_ascii_write_file)
        .def("display", &Histogram_wrapper_display)
        .def("display_survival", &Histogram_wrapper_display_survival)
        .def("spreadsheet_write", &Histogram_wrapper_spreadsheet_write)
        .def("plot_write", &Histogram_wrapper_plot_write)
        .def("plot_write_survival", &Histogram_wrapper_plot_write_survival)
        .def( self == self )
        .def( self != self )
        .def(self_ns::str(self))
    ;

    class_< Parametric, bases< Distribution > , boost::noncopyable, Parametric_Wrapper >("_Parametric", init< const Forward& >())
        .def(init< optional< int, int, int, int, double, double > >())
        .def(init< int, int, int, double, double, optional< double > >())
        .def(init< const Distribution&, optional< int > >())
        .def(init< const Distribution&, double >())
//        .def(init< const Parametric&, double >())
//        .def(init< const Parametric&, optional< char, int > >())
        .def("parametric_mean_computation", &Parametric::parametric_mean_computation)
        .def("parametric_variance_computation", &Parametric::parametric_variance_computation)
        .def("parametric_skewness_computation", &Parametric::parametric_skewness_computation)
        .def("parametric_kurtosis_computation", &Parametric::parametric_kurtosis_computation)
        .def("computation", &Parametric::computation, Parametric_computation_overloads_0_2())
        .def("simulation", &Parametric::simulation)
        .def_readonly("ident", &Parametric::ident)
        .def_readonly("inf_bound", &Parametric::inf_bound)
        .def_readonly("sup_bound", &Parametric::sup_bound)
        .def_readonly("parameter", &Parametric::parameter)
        .def_readonly("probability", &Parametric::probability)
        .def(self_ns::str(self))
    ;


/*    class_< Parametric_model, bases< Parametric > >("Parametricm", init< const Histogram& >() )
        .def(init< const Distribution& >())
        .def(init< const Parametric_model&, optional< bool > >()
        .def(init< const Parametric& >())
        .def(init< const Distribution&, const Histogram& >())
        .def(init< const Histogram& >())
        .def(init< const Parametric&,  const Histogram& >())
        .def("extract_data", &Parametric_model::extract_data, return_value_policy< manage_new_object >())
        .def(self_ns::str(self))
    ;*/

    def("ERROR_LENGTH", ERROR_LENGTH_);
    def("I_DEFAULT", I_DEFAULT_);
    def("D_DEFAULT", D_DEFAULT_);
    def("D_INF", D_INF_);
    def("DOUBLE_ERROR", DOUBLE_ERROR_);
    def("SELF_TRANSITION", SELF_TRANSITION_);

    enum_<UniqueInt<4, 0> >("TestDistributions")
        .value("STANDARD_NORMAL", STANDARD_NORMAL)
        .value("CHI2", CHI2)
        .value("FISHER", FISHER)
        .value("STUDENT", STUDENT)
        .export_values()
    ;

    enum_<UniqueInt<5, 0> >("DistributionIdentifier")
        .value("NON_PARAMETRIC", NONPARAMETRIC)
        .value("BINOMIAL",BINOMIAL)
        .value("POISSON",POISSON)
        .value("NEGATIVE_BINOMIAL",NEGATIVE_BINOMIAL)
        .value("UNIFORM",UNIFORM)
        .export_values()
    ;

    enum_<UniqueInt<4, 1> >("RestorationAlgorithm")
        .value("NO_COMPUTATION", 0)
        .value("FORWARD", FORWARD)
        .value("FORWARD_BACKWARD", FORWARD_BACKWARD)
        .value("VITERBI", VITERBI)
        .export_values()
    ;

    enum_<UniqueInt<6, 0> >("VariableType")
        .value("INT_VALUE", INT_VALUE)
        .value("REAL_VALUE", REAL_VALUE)
        .value("STATE", STATE)
        .value("OLD_INT_VALUE", OLD_INT_VALUE)
        .value("NB_INTERNODE", NB_INTERNODE)
        .value("AUXILIARY", AUXILIARY)
        .export_values()
    ;
}

