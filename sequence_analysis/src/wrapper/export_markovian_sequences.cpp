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
#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;

using namespace sequence_analysis;

#define WRAP MarkovianSequencesWrap

class WRAP {

public:

  static Markovian_sequences*
  markovian_sequences1(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    Markovian_sequences *seq = NULL;
    // TODO
    return seq;
  }

  static Markovian_sequences*
  markovian_sequences2(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    Markovian_sequences *seq = NULL;
    // todo
    return seq;
  }

  static Markovian_sequences*
  markovian_sequences3(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    Markovian_sequences *seq = NULL;
    // todo
    return seq;
  }

  static Distribution_data*
  extract(const Markovian_sequences& input, int type, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, Distribution_data, type,
        variable, value);
  }

  static Markovian_sequences*
  merge(const Markovian_sequences& input, const boost::python::list& input_list)
  {

   CREATE_ARRAY(input_list, const Markovian_sequences*, data);
   SIMPLE_METHOD_TEMPLATE_1(input, merge, Markovian_sequences,
        data_size, data.get());
  }

  static Markovian_sequences*
  cluster_step(const Markovian_sequences& input, int variable, int step,
		  int mode = FLOOR)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, cluster, Markovian_sequences,
    		variable, step, mode);
  }

  static Markovian_sequences*
  cluster_limit(const Markovian_sequences& seq,
		  int variable, boost::python::list& limit, bool add_flag)
  {

     Format_error error;

     int nb_limit = len(limit);
     bool is_float = true;
     int *lint = NULL;
     double *ldouble = NULL;
     Markovian_sequences* ret;
     // Test type
     boost::python::extract<int> get_int(limit[0]);
     if (get_int.check())
       {
         is_float = false;
         lint = new int[nb_limit];
       }
     else
       {
         ldouble = new double[nb_limit];
       }

     // Convert list
     for (int i = 0; i < nb_limit; i++)
       {
         if (is_float)
           {
        	 ldouble[i] = boost::python::extract<double> (limit[i]);
           }
         else
         {
        	 lint[i] = boost::python::extract<int> (limit[i]);
         }
       }

     // Call correct function
     if (is_float)
       {
         ret = seq.cluster(error, variable, nb_limit+1, ldouble);
         delete[] ldouble;
       }
     else
       {
         ret = seq.cluster(error, variable, nb_limit+1, lint, add_flag);
         delete[] lint;
       }

     if (!ret)
       sequence_analysis::wrap_util::throw_error(error);

     return ret;
   }

  static Markovian_sequences*
  transcode(const Markovian_sequences& input, int variable,
      boost::python::list& input_list, bool add_flag = false)
  {

    CREATE_ARRAY(input_list, int, data);
    SIMPLE_METHOD_TEMPLATE_1(input, transcode, Markovian_sequences,
    		variable, data.get(), add_flag);
  }

  static Markovian_sequences*
  merge_variable(const Markovian_sequences& input,
      const boost::python::list& input_seqs,
      int ref_sample)
  {
    CREATE_ARRAY(input_seqs, const Markovian_sequences*, sequences);
    SIMPLE_METHOD_TEMPLATE_1(input, merge_variable, Markovian_sequences,
        sequences_size, sequences.get(), ref_sample);
  }

  static Markovian_sequences*
  select_variable(const Markovian_sequences& input,
      boost::python::list& input_list , bool keep)
  {
    CREATE_ARRAY(input_list, int, data);
    SIMPLE_METHOD_TEMPLATE_1(input, select_variable, Markovian_sequences,
        data_size, data.get(), keep);
  }

  static Markovian_sequences*
  split(const Markovian_sequences& input, int index)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, split, Markovian_sequences, index);
  }

  static Markovian_sequences*
  consecutive_values(const Markovian_sequences& input,
      int variable , bool add_flag )
  {
    HEADER_OS(Markovian_sequences);
    ret = input.consecutive_values(error, os, variable, add_flag);
    FOOTER_OS;
  }

  static Markovian_sequences*
  remove_index_parameter(const Markovian_sequences& input)
  {
     SIMPLE_METHOD_TEMPLATE_0(input, remove_index_parameter, Markovian_sequences);
  }

  static Markovian_sequences*
  add_absorbing_run(const Markovian_sequences& input ,int sequence_length
		, int run_length)
  {
     SIMPLE_METHOD_TEMPLATE_1(input, add_absorbing_run,
    		 Markovian_sequences, sequence_length, run_length);
  }

  static bool
  transition_count(const Markovian_sequences &input, int max_order,
		  bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,
		  const char *path = 0)
  {
    Format_error error;
    std::stringstream os;
    bool ret = true;

    ret = input.transition_count(error, os,
			max_order, begin, estimator);
    FOOTER_OS;
  }

  static std::string
  word_count(const Markovian_sequences& input, int variable, int word_length,
		  int begin_state, int end_state, int min_frequency)
  {
    Format_error error;
    std::stringstream os;
    bool ret = true;

    ret = input.word_count(error, os ,variable, word_length,
        begin_state, end_state,  min_frequency);

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);
    //todo if this what we want to return ?
    return os.str();
  }


  static Variable_order_markov*
  variable_order_markov_estimation1(const Markovian_sequences& input,
      char model_type, int min_order, int max_order, int algorithm,
      double threshold, int estimator, bool global_initial_transition,
      bool global_sample, bool counting_flag)
  {
    HEADER_OS(Variable_order_markov);
    ret = input.variable_order_markov_estimation(error, os, model_type,
        min_order, max_order, algorithm, threshold, estimator,
        global_initial_transition, global_sample, counting_flag);
    FOOTER_OS;
  }

  static Variable_order_markov*
  variable_order_markov_estimation2(const Markovian_sequences& input,
      char type, int max_order, bool global_initial_transition,
      bool counting_flag)
  {
    HEADER(Variable_order_markov);
    ret = input.variable_order_markov_estimation(error, type, max_order,
        global_initial_transition, counting_flag);
    FOOTER;

  }

  static Variable_order_markov*
  variable_order_markov_estimation3(const Markovian_sequences& input,
      const Variable_order_markov &markov, bool global_initial_transition,
      bool counting_flag)
  {
    HEADER_OS(Variable_order_markov);
    ret = input.variable_order_markov_estimation(error, markov,
        global_initial_transition, counting_flag);
    FOOTER_OS;
  }

  static Variable_order_markov*
  lumpability_estimation(const Markovian_sequences& input,
      boost::python::list& input_symbol, int penalty_type, int order,
      bool counting_flag)
  {
    HEADER_OS(Variable_order_markov);

//    int nb_value = len(input_symbol);
//    int *symbol;
  //  symbol = new int[nb_value];

    /*for (int i = 0; i < nb_value; i++)
      {
        symbol[i] = boost::python::extract<int>(input_symbol[i]);
      }
*/
    CREATE_ARRAY(input_symbol, int, symbol);
    ret = input.lumpability_estimation(error, os, symbol.get(), penalty_type,
        order, counting_flag);
    FOOTER_OS;

  }

  static Hidden_variable_order_markov*
  hidden_variable_order_markov_estimation(const Markovian_sequences& input,
      const Hidden_variable_order_markov &hvom, bool global_initial_transition,
      bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(Hidden_variable_order_markov);
    ret = input.hidden_variable_order_markov_estimation(error, os, hvom,
        global_initial_transition, counting_flag, state_sequence, nb_iter);
    FOOTER_OS;
  }


  static Hidden_variable_order_markov*
  hidden_variable_order_markov_stochastic_estimation(
      const Markovian_sequences &input,
      const Hidden_variable_order_markov &hvom, bool global_initial_transition,
      int min_nb_state_sequence, int max_nb_state_sequence, double parameter,
      bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(Hidden_variable_order_markov);
    ret = input.hidden_variable_order_markov_stochastic_estimation(error, os,
        hvom, global_initial_transition, min_nb_state_sequence,
        max_nb_state_sequence, parameter, counting_flag, state_sequence,
        nb_iter);
    FOOTER_OS;
  }

  static bool
  comparison_variable_order_markov(const Markovian_sequences &input,
      boost::python::list &input_markov, char *filename)
  {
    Format_error error;
    bool ret = false;
    std::stringstream os;

    CREATE_ARRAY(input_markov, const Variable_order_markov *, data);
    ret = input.comparison(error, os, data_size, data.get(), filename);
    FOOTER_OS;
  }

  static bool
  comparison_semi_markov(const Markovian_sequences &input,
      boost::python::list &input_markov, char *filename)
  {
    Format_error error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const Semi_markov *, data);
    ret = input.comparison(error, os, data_size, data.get(), filename);
    FOOTER_OS;
  }

  static bool
  comparison_hidden_variable_order_markov(const Markovian_sequences &input,
      boost::python::list &input_markov, int algorithm, const char *filename)
  {
    Format_error error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const Hidden_variable_order_markov *, data);
    ret = input.comparison(error, os, data_size, data.get(), algorithm,
        filename);
    FOOTER_OS;
  }

  static bool
  comparison_hidden_semi_markov(const Markovian_sequences &input,
      boost::python::list &input_markov, int algorithm,
      const char *filename)
  {
    Format_error error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const Hidden_semi_markov*, markov);
    ret = input.comparison(error, os, markov_size, markov.get(), algorithm,
        filename);

    FOOTER_OS;
  }

  static void
  self_transition_computation(Markovian_sequences& input)
  {
      //todo this functio is protected so only the case with no arguments is wrapped for now
      input.self_transition_computation();
  }

  static Hidden_semi_markov*
  hidden_semi_markov_estimation(const Markovian_sequences &input,
      const Hidden_semi_markov &ihsmarkov, int estimator, bool counting_flag,
      bool state_sequence, int nb_iter, int mean_computation)
  {
    HEADER_OS(Hidden_semi_markov);
    ret = input.hidden_semi_markov_estimation(error, os, ihsmarkov, estimator,
        counting_flag, state_sequence, nb_iter, mean_computation);
    FOOTER_OS;
  }

  static Hidden_semi_markov*
  hidden_semi_markov_estimation_model(const Markovian_sequences &input,
      char model_type, int nb_state, bool left_right, int estimator,
      bool counting_flag, bool state_sequence, double occupancy_mean,
      int nb_iter, int mean_computation)
  {
    HEADER_OS(Hidden_semi_markov);
    ret = input.hidden_semi_markov_estimation(error, os, model_type, nb_state,
        left_right, estimator, counting_flag, state_sequence, occupancy_mean,
        nb_iter, mean_computation);
    FOOTER_OS;
  }


  static Hidden_semi_markov*
  hidden_semi_markov_stochastic_estimation(const Markovian_sequences &input,
      const Hidden_semi_markov &ihsmarkov, int min_nb_state_sequence,
      int max_nb_state_sequence, double parameter, int estimator,
      bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(Hidden_semi_markov);

    ret = input.hidden_semi_markov_stochastic_estimation(error, os, ihsmarkov,
        min_nb_state_sequence, max_nb_state_sequence, parameter, estimator,
        counting_flag, state_sequence, nb_iter);

    FOOTER_OS;
  }

  static Hidden_semi_markov*
  hidden_semi_markov_stochastic_estimation_model(
      const Markovian_sequences &input, char model_type, int nb_state,
      bool left_right, int min_nb_state_sequence, int max_nb_state_sequence,
      double parameter, int estimator, bool counting_flag, bool state_sequence,
      double occupancy_mean, int nb_iter)
  {
    HEADER_OS(Hidden_semi_markov);
    ret = input.hidden_semi_markov_stochastic_estimation(error, os, model_type,
        nb_state, left_right, min_nb_state_sequence, max_nb_state_sequence,
        parameter, estimator, counting_flag, state_sequence, occupancy_mean,
        nb_iter);
    FOOTER_OS;

  }

  static Semi_markov*
  semi_markov_estimation(const Markovian_sequences & input, char model_type,
      int estimator, bool counting_flag, int nb_iter, int mean_computation)
  {
    HEADER_OS(Semi_markov);
    ret = input.semi_markov_estimation(error, os, model_type, estimator,
        counting_flag, nb_iter, mean_computation);
    FOOTER_OS;
  }

  static Nonhomogeneous_markov*
  nonhomogeneous_markov_estimation(const Markovian_sequences & input,
      boost::python::list list_ident, bool counting_flag)
  {
    HEADER(Nonhomogeneous_markov);
    CREATE_ARRAY(list_ident, int, data);
    ret = input.nonhomogeneous_markov_estimation(error, data.get(), counting_flag);
    FOOTER;
  }

  static bool
  lumpability_test(const Markovian_sequences & input,
      boost::python::list input_values, int order)
  {
    Format_error error;
    bool ret;
    std::stringstream os;

    CREATE_ARRAY(input_values, int, data);
    ret = input.lumpability_test(error, os, data.get(), order);
    FOOTER_OS;
  }
};

// Boost declaration

void class_markovian_sequences() {




  class_<Markovian_sequences, bases<Sequences> > ("_Markovian_sequences", "Markovian_sequences")

     //todo those constructors if needed
    .def("__init__", make_constructor(WRAP::markovian_sequences1))
    .def("__init__", make_constructor(WRAP::markovian_sequences2))
    .def("__init__", make_constructor(WRAP::markovian_sequences3))

    .def(init <Sequences>())
    .def(init <Markovian_sequences, optional<char, int> >())

    .def(self_ns::str(self)) //__str__

    DEF_RETURN_VALUE("add_absorbing_run", WRAP::add_absorbing_run, args("sequence_length", "run_length"), "addition of a run of absorbing vectors at the end of sequences")
    DEF_RETURN_VALUE("cluster_step", WRAP::cluster_step, args("variable", "step", "mode"), "Cluster")
    DEF_RETURN_VALUE("cluster_limit", WRAP::cluster_limit, args("variable", "limit", "bool_add_flag"), "Cluster")
    DEF_RETURN_VALUE("consecutive_values", WRAP::consecutive_values, args("variable", "AddVariable_flag"), "Consecutive values")
    DEF_RETURN_VALUE("extract", WRAP::extract, args("type", "variable","value"), "Extract distribution data")
    DEF_RETURN_VALUE("split", WRAP::split, args("step"),  "Split")
    DEF_RETURN_VALUE("merge", WRAP::merge, args("sequences"),  "Merge")
    DEF_RETURN_VALUE_NO_ARGS("merge_variable", WRAP::merge_variable, "Merge variables")
    DEF_RETURN_VALUE_NO_ARGS("remove_index_parameter", WRAP::remove_index_parameter, "Remove index parameter")
    DEF_RETURN_VALUE("select_variable", WRAP::select_variable, args("nb_variable", "list variables", "keep"), "select variable")
    DEF_RETURN_VALUE("transcode", WRAP::transcode, args("variable", "symbol", "add_flag"), "Transcode")

    DEF_RETURN_VALUE("hidden_variable_order_markov_estimation", WRAP::hidden_variable_order_markov_estimation, args(""), "todo")
    DEF_RETURN_VALUE("hidden_variable_order_markov_stochastic_estimation", WRAP::hidden_variable_order_markov_stochastic_estimation, args(""), "todo")

    //comparison
    .def("comparison_variable_order_markov", WRAP::comparison_variable_order_markov, args("markov list","filename"), "todo")
    .def("comparison_semi_markov", WRAP::comparison_semi_markov, args("markov list","filename"), "todo")
    .def("comparison_hidden_variable_order_markov", WRAP::comparison_hidden_variable_order_markov, args("markov list","algo","filename"), "todo")
    .def("comparison_hidden_semi_markov", WRAP::comparison_hidden_semi_markov, args("markov list","algo","filename"), "todo")

    //others
    .def("lumpability_test", WRAP::lumpability_test, args("symbols", "order"), "test lumpability and returns status. See LumpabilityTest in sequence_analysis.")
    .def("self_transition_computation", WRAP::self_transition_computation )
    .def("transition_count", WRAP::transition_count, "transition count")
    .def("word_count", WRAP::word_count, args("variable", "word_length","begin_state", "end_state","min_frequency"), "todo" )

    //estimation
    DEF_RETURN_VALUE("variable_order_markov_estimation1", WRAP::variable_order_markov_estimation1, args("model_type", "min_order", "max_order", "algorithm", "threshold", "estimator","global_initial_transition","global_sample", "counting_flag"), "todo")
    DEF_RETURN_VALUE("variable_order_markov_estimation2", WRAP::variable_order_markov_estimation2, args("type","max_order","global_initial_transition","counting_flag"), "todo")
    DEF_RETURN_VALUE("variable_order_markov_estimation3", WRAP::variable_order_markov_estimation3, args("markov","global_initial_transition","counting_flag"), "todo")
    DEF_RETURN_VALUE("lumpability_estimation", WRAP::lumpability_estimation, args("input_symbol", "penalty_type","order","counting_flag"), "todo")
    DEF_RETURN_VALUE("hidden_semi_markov_stochastic_estimation_model", WRAP::hidden_semi_markov_stochastic_estimation_model, args("tobedone"), "todo")
    DEF_RETURN_VALUE("hidden_semi_markov_stochastic_estimation", WRAP::hidden_semi_markov_stochastic_estimation, args("tobedone"), "todo")
    DEF_RETURN_VALUE("hidden_semi_markov_estimation_model", WRAP::hidden_semi_markov_estimation_model, args("tobedone"), "todo")
    DEF_RETURN_VALUE("hidden_semi_markov_estimation", WRAP::hidden_semi_markov_estimation, args("tobedone"), "todo")
    DEF_RETURN_VALUE("semi_markov_estimation", WRAP::semi_markov_estimation, args("tobedone"), "todo")
    DEF_RETURN_VALUE("nonhomogeneous_markov_estimation", WRAP::nonhomogeneous_markov_estimation, args("tobedone"), "todo")
    ;


  /*
   Markovian_sequences* remove_variable_1() const;
   Markovian_sequences* initial_run_computation(Format_error &error) const;

   std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,   bool exhaustive = false) const;
   bool ascii_data_write(Format_error &error , const char *path , char format = 'c' , bool exhaustive = false) const;
   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(Format_error &error , const char *path ,   bool exhaustive = false) const;
   bool spreadsheet_write(Format_error &error , const char *path) const;
   bool plot_write(Format_error &error , const char *prefix ,   const char *title = 0) const;

   bool mtg_write(Format_error &error , const char *path , int *itype) const;

   double iid_information_computation() const;

   void sojourn_time_histogram_computation(int variable);
   void create_observation_histogram(int nb_state);
   void observation_histogram_computation();
   void build_observation_histogram();
   void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,  bool initial_run_flag = false);

   // acces membres de la classe

   Curves* get_self_transition(int state) const { return self_transition[state]; }
   Histogram*** get_observation() const { return observation; }
   Histogram** get_observation(int variable) const { return observation[variable]; }
   Histogram* get_observation(int variable , int state) const   { return observation[variable][state]; }
   Sequence_characteristics* get_characteristics(int variable) const   { return characteristics[variable]; }
   */

}

#undef WRAP


void class_self_transition() {

  class_<Self_transition, bases<Curves> > ("_Self_transition", "Self_transition", no_init)
    .def(init <int>())
    DEF_RETURN_VALUE_NO_ARGS("monomolecular_regression", &Self_transition::monomolecular_regression, "returns monomolecular regression")
    DEF_RETURN_VALUE_NO_ARGS("logistic_regression", &Self_transition::logistic_regression,"returns logistic regression")
    ;
    //Done
}

