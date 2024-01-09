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
 *        $Id: export_distribution.cpp 18034 2015-04-23 07:15:31Z guedon $
 *
 *-----------------------------------------------------------------------------*/



#include "wrapper_util.h"
#include "export_distribution.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"


#include <boost/python.hpp>
// definition of boost::python::len
#include <boost/python/detail/api_placeholder.hpp>
// definition of boost::python::make_constructor
#include <boost/python/make_constructor.hpp>
#include <boost/shared_ptr.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;


#include "boost_python_aliases.h"

////////////////////// Export DiscreteParametricModel ////////////////////////////////



class DistributionWrap {

public:

  static MultiPlotSet*
  get_plotable_dists(const Distribution &p,
      const boost::python::list& dist_list)
  {
    StatError error;
    int nb_dist = boost::python::len(dist_list);
    stat_tool::wrap_util::auto_ptr_array<const Distribution *> dists(
        new const Distribution*[nb_dist]);

    for (int i = 0; i < nb_dist; i++)
      dists[i] = extract<const Distribution*> (dist_list[i]);

    const Distribution** d = dists.get();

    MultiPlotSet *ret = p.get_plotable_distributions(error, nb_dist, d);
    if (!ret)
      stat_tool::wrap_util::throw_error(error);
//    for (int i = 0; i < nb_dist; i++)
 //       delete dists[i];

    return ret;
  }

  static MultiPlotSet*
  survival_get_plotable(const Distribution &p)
  {
    StatError error;
    MultiPlotSet *ret = p.survival_get_plotable(error);
    if (!ret)
      ERROR;
    return ret;
  }

  static MultiPlotSet*
  get_plotable(const Distribution &p)
  {
    MultiPlotSet *ret = p.get_plotable();
    return ret;
  }

    // survival_ascii_write wrapping
    WRAP_METHOD_SURVIVAL_ASCII_WRITE( Distribution);

    //survival_spreadsheet_write wrapping
    WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE( Distribution);

    // survival_plot_write wrapping
    WRAP_METHOD_SURVIVAL_PLOT_WRITE( Distribution);

    //truncate
    WRAP_METHOD1(Distribution, truncate, DiscreteParametricModel, int);

};




// Boost.Python Wrapper export function
void class_distribution()
{

#define WRAP DistributionWrap
  // Distribution base class
  class_< Distribution>("_Distribution")
    .def(init< boost::python::optional< int > > ())
    .def(init<const FrequencyDistribution&>())
    .def(init<const Distribution&, double>())
    .def(init<const Distribution&, boost::python::optional< char, int > >())

    .def(self_ns::str(self)) // __str__
    .def( self == self )
    .def( self != self )


    // todo siwth cto properties ?
    .add_property("nb_value", &Distribution::nb_value,
        "Number of values above zero")
    .def_readonly("get_alloc_nb_value", &Distribution::alloc_nb_value,
        "Number of values with zero probability")
    .def_readonly("get_max", &Distribution::max,
        "probability maximum")
    .def_readonly("get_complement", &Distribution::complement,
        "complementary probability")
    .def_readonly("get_mean", &Distribution::mean,
        "mean")
    .def_readonly("get_variance", &Distribution::variance,
        "variance")
    .def_readonly("get_nb_parameter", &Distribution::nb_parameter,
        "number of unknown parameters")

    .def("simulation", &Distribution::simulation, "simulate one realization")

    // no tested. is it useful ?
    DEF_RETURN_VALUE_NO_ARGS("get_plotable_list", WRAP::get_plotable_dists,
        "Return a plotable for a list of distribution")
    DEF_RETURN_VALUE_NO_ARGS("survival_get_plotable", WRAP::survival_get_plotable,
        "Return a survival plotable")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable,
        "Return a plotable")

    .def("survival_ascii_write", WRAP::survival_ascii_write,
        "Return a string containing the object description (survival viewpoint)")
    .def("survival_plot_write", WRAP::survival_plot_write,
        args("prefix", "title"), "Write GNUPLOT files (survival viewpoint)")
    .def("survival_spreadsheet_write", WRAP::survival_spreadsheet_write,
        args("filename"),"Write object to filename (spreadsheet format)")
    DEF_RETURN_VALUE("truncate", WRAP::truncate,
        args("index"),  "See Truncate")


    /*
     *
   NOT DONE but may be considered::

   double *mass;           // probabilites de chaque valeur
   double *cumul;          // fonction de repartition

   void max_computation();
   void mean_computation();
   void variance_computation();

   std::ostream& ascii_characteristic_print(...
   std::ostream& spreadsheet_characteristic_print(...
   std::ostream& spreadsheet_print(...
   std::ostream& spreadsheet_print(

   int plot_nb_value_computation(const FrequencyDistribution *histo = 0) const;

   bool plot_print(const char *path , double *concentration , double scale) const;
   bool plot_print(const char *path , const FrequencyDistribution *histo = 0) const;

   virtual std::ostream& plot_title_print(std::ostream &os) const { return os; }
   bool survival_plot_print(const char *path , double *survivor) const;
   std::ostream& print(std::ostream&) const;

   void plotable_mass_write(SinglePlot &plot , double scale = 1.) const;
   void plotable_cumul_write(SinglePlot &plot) const;
   void plotable_cumul_matching_write(SinglePlot &plot , const Distribution &reference_dist) const;
   void plotable_concentration_write(SinglePlot &plot) const;
   void plotable_survivor_write(SinglePlot &plot) const;

   void convolution(Distribution &dist1 , Distribution &dist2 ,  int inb_value = I_DEFAULT);

   void nb_value_computation();
   void offset_computation();
   double concentration_computation() const;

   void cumul_computation();
   double* survivor_function_computation() const;
   double* concentration_function_computation() const;
   void log_computation();

   double survivor_likelihood_computation(const FrequencyDistribution &histo) const;
   double chi2_value_computation(const FrequencyDistribution &histo) const;
   void chi2_degree_of_freedom(const FrequencyDistribution &histo , Test &test) const;

   void penalty_computation(double weight , int type , double *penalty , int outside) const;


   bool plot_write(StatError &error , const char *prefix , int nb_dist ,
                   const Distribution **idist , const char *title) const;



   double mean_absolute_deviation_computation() const;
   double skewness_computation() const;
   double kurtosis_computation() const;
   double information_computation() const;

   double first_difference_norm_computation() const;
   double second_difference_norm_computation() const;

   double likelihood_computation(const Reestimation<int> &histo) const   { return histo.likelihood_computation(*this); }
   double likelihood_computation(const Reestimation<double> &histo) const   { return histo.likelihood_computation(*this); }
   void chi2_fit(const FrequencyDistribution &histo , Test &test) const;


   DiscreteParametricModel* truncate(StatError &error , int imax_value) const;
*/


    ;
#undef WRAP

}




class DiscreteParametricWrap
{

public:



};

void class_discrete_parametric()
{

#define WRAP DiscreteParametricWrap

  // DiscreteParametric base class
  class_< DiscreteParametric, bases< Distribution > >
    ("_DiscreteParametric", init< boost::python::optional< int, int, int, int, double, double > >())
    .def(init<int, int, int, double, double, boost::python::optional< double > >())
    .def(init<int, int>())
    .def(init<const Distribution&, boost::python::optional<int> >())
    .def(init<const Distribution&, double>())
    .def(init<const DiscreteParametric&, double>())
    .def(init<const FrequencyDistribution& >())
    .def(init<const DiscreteParametric&, boost::python::optional< char, int> >())
    .def(init<Distribution&>())
    .def(init<DiscreteParametric&>())
    //to remove
    .def(init<Distribution&>())
    .def_readonly("get_ident", &DiscreteParametric::ident)
    .def_readonly("get_inf_bound", &DiscreteParametric::inf_bound)
    .def_readonly("get_sup_bound", &DiscreteParametric::sup_bound)
    .def_readonly("get_parameter", &DiscreteParametric::parameter)
    .def_readonly("get_probability", &DiscreteParametric::probability)
    .def(self_ns::str(self))
    .def("simulate", &DiscreteParametric::simulation, "Simulation one value")

    ;

  /*
   *
   NOT DONE but may be considered later

    std::ostream& ascii_print(std::ostream &os) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;
    std::ostream& plot_title_print(std::ostream &os) const;

    void nb_parameter_update();

    void binomial_computation(int inb_value , char mode);
    void poisson_computation(int inb_value , double cumul_threshold , char mode);
    void negative_binomial_computation(int inb_value , double cumul_threshold , char mode);
    void uniform_computation();

    double renewal_likelihood_computation(...
    void expectation_step(...
    void reestimation(const Reestimation<double> *reestim , int nb_estim = 1);

    double state_occupancy_likelihood_computation(...
    double state_occupancy_likelihood_computation(...
    void expectation_step(...
    void expectation_step(...
    int nb_parameter_computation();

    double parametric_mean_computation() const;
    double parametric_variance_computation() const;
    double parametric_skewness_computation() const;
    double parametric_kurtosis_computation() const;

    void computation(int min_nb_value = 1 , double cumul_threshold = CUMUL_THRESHOLD);
  */

#undef WRAP

}





// Wrapper class

class DiscreteParametricModelWrap
{

public:

  static boost::shared_ptr<DiscreteParametricModel> parametric_model_from_file(char* filename)
  {
    StatError error;
    DiscreteParametricModel *model = NULL;
    model = discrete_parametric_ascii_read(error, filename);

    if(!model) stat_tool::wrap_util::throw_error(error);
    return boost::shared_ptr<DiscreteParametricModel>(model);
  }


  // simulation method wrapping
  WRAP_METHOD1(DiscreteParametricModel, simulation, DiscreteDistributionData, int);

  // extract_data method wrapping
  WRAP_METHOD0(DiscreteParametricModel, extract_data, DiscreteDistributionData);

  // survival_ascii_write wrapping
  WRAP_METHOD_SURVIVAL_ASCII_WRITE(DiscreteParametricModel);

  //survival_spreadsheet_write wrapping
  WRAP_METHOD_SURVIVAL_SPREADSHEET_WRITE(DiscreteParametricModel);

  // survival_plot_write wrapping
  WRAP_METHOD_SURVIVAL_PLOT_WRITE(DiscreteParametricModel);

//   static void plot_write(const DiscreteParametricModel& p,
// 			 const std::string& prefix, const std::string& title,
// 			 const boost::python::list& dist_list)
//   {
//     StatError error;

//     int nb_dist = boost::python::len(dist_list);
//     stat_tool::wrap_util::auto_ptr_array<const Distribution *>
//       dists(new const Distribution*[nb_dist]);

//     const Distribution &d = (const Distribution&)(p);

//     if(!d.plot_write(error, prefix.c_str(), nb_dist, dists.get(), title.c_str()))
//       stat_tool::wrap_util::throw_error(error);
//   }


  static MultiPlotSet* get_plotable(const DiscreteParametricModel& p,
				const boost::python::list& dist_list)
  {

    cout << "get_plotable" << endl;
    StatError error;
    int nb_dist = boost::python::len(dist_list);
    stat_tool::wrap_util::auto_ptr_array<const Distribution *>
      dists(new const Distribution*[nb_dist]);

    for (int i = 0; i < nb_dist; i++)
      dists[i] = extract<const Distribution*>(dist_list[i]);

    const Distribution** d = dists.get();

    MultiPlotSet* ret = p.get_plotable_distributions(error, nb_dist, d);
    if(!ret)
      stat_tool::wrap_util::throw_error(error);

    //for (int i = 0; i < nb_dist; i++)
    //    delete dists[i];

    return ret;
  }

  //survival_get_plotable wrapping
  WRAP_METHOD_SURVIVAL_GET_PLOTABLE(DiscreteParametricModel);

  //file_ascii_write wrapping
  WRAP_METHOD_FILE_ASCII_WRITE(DiscreteParametricModel);

};




void class_discrete_parametric_model()
{

#define WRAP DiscreteParametricModelWrap


  // DiscreteParametricModel
  class_< DiscreteParametricModel, bases< DiscreteParametric, StatInterface > >
  ("_DiscreteParametricModel", "Parametric model", init <const FrequencyDistribution& >())

  .def(init< int, int, int, double, double, boost::python::optional< double > >())
  // this constructor clashes with the previous one and fail to pass tests.
  //.def(init< int, boost::python::optional <int, int, int, double, double > >())
  .def(init <const Distribution& >())
  .def(init <const DiscreteParametric& >())
  .def(init <const DiscreteParametricModel& , boost::python::optional< bool> >())

  .def("__init__", make_constructor(DiscreteParametricModelWrap::parametric_model_from_file))
  .def(self_ns::str(self)) // __str__

  /*.def("get_histogram", &DiscreteParametricModel::get_histogram,
    	return_value_policy< manage_new_object >(),
    	"returns histogram")
*/
    // Output
    .def("get_plotable", DiscreteParametricModelWrap::get_plotable,
    		return_value_policy< manage_new_object >(),
    		"Return a plotable for a list of distribution")
    .def("get_plotable", &StatInterface::get_plotable,
    		return_value_policy< manage_new_object >(),
    		"Return a plotable (no parameters)")

//     .def("plot_write", DiscreteParametricModelWrap::plot_write,
// 	 args("prefix", "title", "dists"),
// 	 "Write GNUPLOT files (with prefix) for a list of distribution")

//     .def("plot_write", &StatInterfaceWrap::plot_write,
// 	 args("prefix", "title"),
// 	  "Write GNUPLOT files (with prefix)")

    .def("survival_ascii_write", WRAP::survival_ascii_write,
    		"Return a string containing the object description (survival viewpoint)")
    .def("survival_plot_write", WRAP::survival_plot_write,
    		args("prefix", "title"),
    		"Write GNUPLOT files (survival viewpoint)")
    .def("survival_get_plotable", WRAP::survival_get_plotable,
    		return_value_policy< manage_new_object >(),
    		"Return a plotable object")
    .def("survival_spreadsheet_write", WRAP::survival_spreadsheet_write,
    		args("filename"),
    		"Write object to filename (spreadsheet format)")
    .def("extract_data", WRAP::extract_data,
    		return_value_policy< manage_new_object >(),
    		"Return the 'data' part of the model")
    .def("simulate", WRAP::simulation,
    		return_value_policy< manage_new_object >(),
    		args("nb_value"),
    		"Simulate values")
    .def("simulate", &DiscreteParametric::simulation,
    		"Simulate one value")
    .def("file_ascii_write", WRAP::file_ascii_write,
    		"Return a string containing the object description")
    ;

  //remains to be done if needed
  /*
    DiscreteParametricModel(const Distribution &dist , const FrequencyDistribution *histo);
    DiscreteParametricModel(const DiscreteParametric &dist , const FrequencyDistribution *histo);

    std::ostream& line_write(std::ostream &os) const;
    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,const char *title = 0) const;

*/
#undef WRAP
}
