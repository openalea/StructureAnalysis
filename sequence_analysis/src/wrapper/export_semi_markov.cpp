/*------------------------------------------------------------------------------
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guédon <yann.guedon@cirad.fr>
 *                        Thomas Cokelaer <Thomas.Cokelaer@inria.fr>
 *
 *        Distributed under the GPL 2.0 License.
 *        See accompanying file LICENSE.txt or copy at
 *           http://www.gnu.org/licenses/gpl-2.0.txt
 *
 *        OpenAlea WebSite : http://openalea.gforge.inria.fr
 *
 *        $Id:  $
 *
 *-----------------------------------------------------------------------------*/

#include "wrapper_util.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/sequence_label.h"
#include "tool/config.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>
#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;

#define WRAP SemiMarkovWrap
class SemiMarkovWrap
{

public:

  static boost::shared_ptr<Semi_markov>
  semi_markov_from_file(char* filename, int length, bool counting_flag, double cumul_threshold)
  {
    Format_error error;

    Semi_markov *semi_markov = NULL;

    semi_markov = semi_markov_ascii_read(error, filename, length,
        counting_flag, cumul_threshold);


    return boost::shared_ptr<Semi_markov>(semi_markov);
  }

  static Semi_markov_data*
  extract_data(const Semi_markov& input)
  {
    Format_error err;
    Semi_markov_data *ret;
    ret = input.extract_data(err);
    return ret;
    //SIMPLE_METHOD_TEMPLATE_1(input, extract_data, Semi_markov_data);
  }

  static Parametric_model*
  extract_histogram(const Semi_markov& input, int state, int histogram_type)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, Parametric_model, state, histogram_type);
  }

  static Parametric_model*
  extract(const Semi_markov& input, int type, int variable, int value)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, extract, Parametric_model, type, variable,
        value);
  }



  static Parametric_model*
  get_forward(const Semi_markov& input)
  {
    // todo: output must be list of Forward
    // todo: export Forward in stat_tool!
    Parametric_model *ret = NULL;
    return ret;
  }

  static Forward*
  get_forward(const Semi_markov& input, int state)
  {
    Forward *ret;
    ret = input.get_forward(state);
    return ret;
  }
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_forward_overloads, SemiMarkovWrap::get_forward, 0, 1);


  static void
  file_ascii_write(const Semi_markov& d, const char* path, bool exhaustive)
  {
    bool result = true;
    Format_error error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }



  static Semi_markov*
  thresholding(const Semi_markov& input, double min_probability)
  {
	return input.thresholding(min_probability);
  }

  static Semi_markov_data*
  simulation_histogram(const Semi_markov& input,
	   const Histogram &hlength , bool counting_flag, bool divergence_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, Semi_markov_data,
   	    hlength, counting_flag, divergence_flag);
  }

  static Semi_markov_data*
  simulation_nb_elements(const Semi_markov& input,
		 int nb_sequence , int length , bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, Semi_markov_data,
   	    nb_sequence, length, counting_flag);

  }

  static Semi_markov_data*
  simulation_markovian_sequences(const Semi_markov& input,
		  int nb_sequence , const Markovian_sequences &iseq , bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, Semi_markov_data,
   	    nb_sequence, iseq, counting_flag);
  }


  static Distance_matrix*
  divergence_computation_histo(const Semi_markov &input,
      boost::python::list input_markov,
      boost::python::list input_histogram_length,  const char *filename)
  {
    HEADER_OS(Distance_matrix);
    CREATE_ARRAY(input_markov, const Semi_markov *, markov);
    CREATE_ARRAY(input_histogram_length, Histogram *, hlength);
    ret = input.divergence_computation(error, os, markov_size, markov.get(), hlength.get(), filename);
    FOOTER_OS;
  }

  static Distance_matrix*
  divergence_computation_length(const Semi_markov &input,
     boost::python::list input_markov , int nb_sequence ,
     int length , const char *filename)
  {
    HEADER_OS(Distance_matrix);
    CREATE_ARRAY(input_markov, const Semi_markov *, markov);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_sequence, length, filename);
    FOOTER_OS;
  }


  static Distance_matrix*
  divergence_computation_sequences(const Semi_markov &input,
      boost::python::list &input_markov, boost::python::list &input_sequences,
      int nb_seq, const char *filename)
  {
    // there is as much Variable_order_markov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    HEADER_OS(Distance_matrix);

    CREATE_ARRAY(input_markov, const Semi_markov *, markov);
    CREATE_ARRAY(input_sequences, const Markovian_sequences *, sequences);

    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_seq, sequences.get(), filename);
    FOOTER_OS;
  }

};

// Boost declaration

void
class_semi_markov()
{

  class_<Semi_markov, bases<STAT_interface> >
  ("_Semi_markov",  "Semi_markov\n"
      "Constructors from a file required 3 arguments: length(int), counting_flag(boolean) and cumul_threshold (double)")
    .def("__init__", make_constructor(SemiMarkovWrap::semi_markov_from_file))

    .def(self_ns::str(self)) //__str__

    .add_property("nb_iterator", &Semi_markov::get_nb_iterator, "returns nb iterator")
    .add_property("nb_output_process", &Semi_markov::get_nb_output_process, "returns nb output process")

    .def("get_state_subtype", &Semi_markov::get_state_subtype, args("index"), "returns state subtype")

    //DEF_RETURN_VALUE("get_parametric_process", &Semi_markov::get_parametric_process, args("variable_index"), "returns parametric process corresponding to the given variable")
    .def("extract_histogram", SemiMarkovWrap::extract_histogram, return_value_policy< manage_new_object >(), "todo")
    .def("extract", SemiMarkovWrap::extract, return_value_policy< manage_new_object >(), "todo")

    .def("get_forward", (Forward *(*)(const Semi_markov&, int)) SemiMarkovWrap::get_forward, return_value_policy< manage_new_object >(), SemiMarkovWrap::get_forward_overloads())

    DEF_RETURN_VALUE_NO_ARGS("get_semi_markov_data", &Semi_markov::get_semi_markov_data, "returns semi_markov_data")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", SemiMarkovWrap::extract_data, "returns semi_markov_data")

    .def("file_ascii_write", SemiMarkovWrap::file_ascii_write,"Save vector summary into a file")
    DEF_RETURN_VALUE("thresholding", SemiMarkovWrap::thresholding, args("index"), "todo")

    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_nb_elements",WRAP::simulation_nb_elements, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("todo"), "simulation")

    DEF_RETURN_VALUE("divergence_computation_histo", WRAP::divergence_computation_histo, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_length", WRAP::divergence_computation_length, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_sequences", WRAP::divergence_computation_sequences, args("input", "input_markov", "input_sequence", "filename"), "todo")

    ;


  //todo file_ascii_write not accessible ?
  /*
   *
   Semi_markov();
   Semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value);
   Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy , const Nonparametric_process *pobservation , int length ,  bool counting_flag);
   Semi_markov(const Semi_markov &smarkov , bool data_flag = true ,  int param = I_DEFAULT)  :Chain(smarkov) { copy(smarkov , data_flag , param); }

    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,const char *title = 0) const;

   void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
   void characteristic_computation(const Semi_markov_data &seq , bool counting_flag, int variable = I_DEFAULT , bool length_flag = true);

   double likelihood_computation(const Markovian_sequences &seq , int index) const;
   double likelihood_computation(const Semi_markov_data &seq) const;


   Forward** get_forward() const { return forward; }

   Nonparametric_sequence_process* get_nonparametric_process(int variable)   const { return nonparametric_process[variable]; }
   Parametric_process** get_parametric_process() const { return parametric_process; }


   */

  ;
}
#undef WRAP


class SemiMarkovDataWrap
{

public:

   static Semi_markov_data*
   remove_index_parameter(const Semi_markov_data& input)
   {
     SIMPLE_METHOD_TEMPLATE_0(input, remove_index_parameter, Semi_markov_data);
   }

   static Distribution_data*
   extract(const Semi_markov_data& input, int type, int variable, int value)
   {
     SIMPLE_METHOD_TEMPLATE_1(input, extract, Distribution_data, type, variable, value);
   }

};

void
class_semi_markov_data()
{

  class_<Semi_markov_data, bases<Markovian_sequences> > ("_Semi_markov_data", "Semi_markov_data")
    .def(init<Markovian_sequences> ())
    .def(init< Markovian_sequences, char, bool> ())
    .def(init<Semi_markov_data, bool, char> ())

    .add_property("likelihood", &Semi_markov_data::get_likelihood, "returns likelihood")
    .add_property("hidden_likelihood", &Semi_markov_data::get_hidden_likelihood, "returns hidden likelihood")

    .def("get_posterior_probability", &Semi_markov_data::get_posterior_probability, args("index"))

    DEF_RETURN_VALUE_NO_ARGS("get_semi_markov", &Semi_markov_data::get_semi_markov, "returns semi_markov")
    DEF_RETURN_VALUE_NO_ARGS("get_chain_data", &Semi_markov_data::get_chain_data, "returns chain data")
    DEF_RETURN_VALUE_NO_ARGS("remove_index_parameter", &SemiMarkovDataWrap::remove_index_parameter, "remove index parameter")

 ;

  // DONE
  // todo do we need the constructor and build_transition?
  /*
   Semi_markov_data();
   Semi_markov_data(const Histogram &ihlength , int inb_variable , bool init_flag = false);
   void build_transition_count(const Semi_markov *smarkov = 0);
   */

}

