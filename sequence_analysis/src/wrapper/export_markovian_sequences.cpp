/*------------------------------------------------------------------------------
 *
 *        VPlants.Sequence_analysis : VPlants Statistics module
 *
 *        Copyright 2006-2007 INRIA - CIRAD - INRA
 *
 *        File author(s): Yann Guedon <yann.guedon@cirad.fr>
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
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/sequence_label.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

#include "boost_python_aliases.h"

using namespace boost::python;
using namespace boost;
using namespace stat_tool;
using namespace sequence_analysis;



#define WRAP MarkovianSequencesWrap

class WRAP {

public:

  static MarkovianSequences*
  markovian_sequences1(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    MarkovianSequences *seq = NULL;
    // TODO
    return seq;
  }

  static MarkovianSequences*
  markovian_sequences2(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    MarkovianSequences *seq = NULL;
    // todo
    return seq;
  }

  static MarkovianSequences*
  markovian_sequences3(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    MarkovianSequences *seq = NULL;
    // todo
    return seq;
  }

  static DiscreteDistributionData*
  extract(const MarkovianSequences &input, int type, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteDistributionData, type,
        variable, value);
  }

  static MarkovianSequences*
  merge(const MarkovianSequences &input, const boost::python::list& input_list)
  {

   CREATE_ARRAY(input_list, const MarkovianSequences*, data);
   SIMPLE_METHOD_TEMPLATE_1(input, merge, MarkovianSequences,
        data_size, data.get());
  }

  static MarkovianSequences*
  cluster_step(const MarkovianSequences &input, int variable, int step,
		  int mode = FLOOR)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, cluster, MarkovianSequences,
    		variable, step, mode);
  }

  static MarkovianSequences*
  cluster_limit(const MarkovianSequences &seq,
		  int variable, boost::python::list& limit, bool add_flag)
  {

     StatError error;

     int nb_limit = len(limit);
     bool is_float = true;
     int *lint = NULL;
     double *ldouble = NULL;
     MarkovianSequences* ret;
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

  static MarkovianSequences*
  transcode(const MarkovianSequences &input, int variable,
      boost::python::list& input_list, bool add_flag = false)
  {

    CREATE_ARRAY(input_list, int, data);
    SIMPLE_METHOD_TEMPLATE_1(input, transcode, MarkovianSequences,
    		variable, data.get(), add_flag);
  }

  static MarkovianSequences*
  merge_variable(const MarkovianSequences &input,
      const boost::python::list& input_seqs,
      int ref_sample)
  {
    CREATE_ARRAY(input_seqs, const MarkovianSequences*, sequences);
    SIMPLE_METHOD_TEMPLATE_1(input, merge_variable, MarkovianSequences,
        sequences_size, sequences.get(), ref_sample);
  }

  static MarkovianSequences*
  select_variable(const MarkovianSequences &input,
      boost::python::list& input_list , bool keep)
  {
    CREATE_ARRAY(input_list, int, data);
    SIMPLE_METHOD_TEMPLATE_1(input, select_variable, MarkovianSequences,
        data_size, data.get(), keep);
  }

  static MarkovianSequences*
  split(const MarkovianSequences &input, int index)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, split, MarkovianSequences, index);
  }

  static MarkovianSequences*
  consecutive_values(const MarkovianSequences &input,
      int variable , bool add_flag )
  {
    HEADER_OS(MarkovianSequences);
    ret = input.consecutive_values(error, os, variable, add_flag);
    FOOTER_OS;
  }

  static MarkovianSequences*
  remove_index_parameter(const MarkovianSequences &input)
  {
     SIMPLE_METHOD_TEMPLATE_0(input, remove_index_parameter, MarkovianSequences);
  }

  static MarkovianSequences*
  add_absorbing_run(const MarkovianSequences &input ,int sequence_length
		, int run_length)
  {
     SIMPLE_METHOD_TEMPLATE_1(input, add_absorbing_run,
    		 MarkovianSequences, sequence_length, run_length);
  }

  static bool
  transition_count(const MarkovianSequences &input, int max_order,
		  bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,
		  const char *path = 0)
  {
    StatError error;
    std::stringstream os;
    bool ret = true;

    ret = input.transition_count(error, os,
			max_order, begin, estimator);
    FOOTER_OS;
  }

  static std::string
  word_count(const MarkovianSequences &input, int variable, int word_length,
		  int begin_state, int end_state, int min_frequency)
  {
    StatError error;
    std::stringstream os;
    bool ret = true;

    ret = input.word_count(error, os ,variable, word_length,
        begin_state, end_state,  min_frequency);

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);
    //todo if this what we want to return ?
    return os.str();
  }


  static VariableOrderMarkov*
  variable_order_markov_estimation1(const MarkovianSequences &input,
      char model_type, int min_order, int max_order, int algorithm,
      double threshold, int estimator, bool global_initial_transition,
      bool global_sample, bool counting_flag)
  {
    HEADER_OS(VariableOrderMarkov);
    ret = input.variable_order_markov_estimation(error, os, model_type,
        min_order, max_order, algorithm, threshold, estimator,
        global_initial_transition, global_sample, counting_flag);
    FOOTER_OS;
  }

  static VariableOrderMarkov*
  variable_order_markov_estimation2(const MarkovianSequences &input,
      char type, int max_order, bool global_initial_transition,
      bool counting_flag)
  {
    HEADER(VariableOrderMarkov);
    ret = input.variable_order_markov_estimation(error, type, max_order,
        global_initial_transition, counting_flag);
    FOOTER;

  }

  static VariableOrderMarkov*
  variable_order_markov_estimation3(const MarkovianSequences &input,
      const VariableOrderMarkov &markov, bool global_initial_transition,
      bool counting_flag)
  {
    HEADER_OS(VariableOrderMarkov);
    ret = input.variable_order_markov_estimation(error, markov,
        global_initial_transition, counting_flag);
    FOOTER_OS;
  }

  static VariableOrderMarkov*
  lumpability_estimation(const MarkovianSequences &input,
      boost::python::list& input_symbol, int penalty_type, int order,
      bool counting_flag)
  {
    HEADER_OS(VariableOrderMarkov);

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

  static HiddenVariableOrderMarkov*
  hidden_variable_order_markov_estimation(const MarkovianSequences &input,
      const HiddenVariableOrderMarkov &hvom, bool global_initial_transition,
      bool common_dispersion, bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(HiddenVariableOrderMarkov);
    ret = input.hidden_variable_order_markov_estimation(error, os, hvom,
        global_initial_transition, common_dispersion, counting_flag,
        state_sequence, nb_iter);
    FOOTER_OS;
  }


  static HiddenVariableOrderMarkov*
  hidden_variable_order_markov_stochastic_estimation(
      const MarkovianSequences &input,
      const HiddenVariableOrderMarkov &hvom, bool global_initial_transition,
      bool common_dispersion, int min_nb_state_sequence, int max_nb_state_sequence,
      double parameter, bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(HiddenVariableOrderMarkov);
    ret = input.hidden_variable_order_markov_stochastic_estimation(error, os,
        hvom, global_initial_transition, common_dispersion, min_nb_state_sequence,
        max_nb_state_sequence, parameter, counting_flag, state_sequence,
        nb_iter);
    FOOTER_OS;
  }

  static bool
  comparison_variable_order_markov(const MarkovianSequences &input,
      boost::python::list &input_markov, char *filename)
  {
    StatError error;
    bool ret = false;
    std::stringstream os;

    CREATE_ARRAY(input_markov, const VariableOrderMarkov*, data);
    ret = input.comparison(error, os, data_size, data.get(), filename);
    FOOTER_OS;
  }

  static bool
  comparison_semi_markov(const MarkovianSequences &input,
      boost::python::list &input_markov, char *filename)
  {
    StatError error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const SemiMarkov*, data);
    ret = input.comparison(error, os, data_size, data.get(), filename);
    FOOTER_OS;
  }

  static bool
  comparison_hidden_variable_order_markov(const MarkovianSequences &input,
      boost::python::list &input_markov, int algorithm, const char *filename)
  {
    StatError error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const HiddenVariableOrderMarkov*, data);
    ret = input.comparison(error, os, data_size, data.get(), algorithm,
        filename);
    FOOTER_OS;
  }

  static bool
  comparison_hidden_semi_markov(const MarkovianSequences &input,
      boost::python::list &input_markov, int algorithm,
      const char *filename)
  {
    StatError error;
    std::stringstream os;
    bool ret = false;

    CREATE_ARRAY(input_markov, const HiddenSemiMarkov*, markov);
    ret = input.comparison(error, os, markov_size, markov.get(), algorithm,
        filename);

    FOOTER_OS;
  }

  static void
  self_transition_computation(MarkovianSequences &input)
  {
      //todo this functio is protected so only the case with no arguments is wrapped for now
      input.self_transition_computation();
  }

  static HiddenSemiMarkov*
  hidden_semi_markov_estimation(const MarkovianSequences &input,
      const HiddenSemiMarkov &ihsmarkov, bool common_dispersion, int estimator,
      bool counting_flag, bool state_sequence, int nb_iter, int mean_computation_method)
  {
    HEADER_OS(HiddenSemiMarkov);
    ret = input.hidden_semi_markov_estimation(error, os, ihsmarkov, common_dispersion,
        estimator, counting_flag, state_sequence, nb_iter, mean_computation_method);
    FOOTER_OS;
  }

  static HiddenSemiMarkov*
  hidden_semi_markov_estimation_model(const MarkovianSequences &input,
      char model_type, int nb_state, bool left_right, double occupancy_mean,
      bool common_dispersion, int estimator, bool counting_flag,
      bool state_sequence, int nb_iter, int mean_computation_method, bool random_initialization)
  {
    HEADER_OS(HiddenSemiMarkov);
    ret = input.hidden_semi_markov_estimation(error, os, model_type, nb_state,
        left_right, occupancy_mean, common_dispersion, estimator, counting_flag,
        state_sequence, nb_iter, mean_computation_method, random_initialization);
    FOOTER_OS;
  }


  static HiddenSemiMarkov*
  hidden_semi_markov_stochastic_estimation(const MarkovianSequences &input,
      const HiddenSemiMarkov &ihsmarkov, bool common_dispersion, int min_nb_state_sequence,
      int max_nb_state_sequence, double parameter, int estimator,
      bool counting_flag, bool state_sequence, int nb_iter)
  {
    HEADER_OS(HiddenSemiMarkov);

    ret = input.hidden_semi_markov_stochastic_estimation(error, os, ihsmarkov,
        common_dispersion, min_nb_state_sequence, max_nb_state_sequence, parameter,
        estimator, counting_flag, state_sequence, nb_iter);

    FOOTER_OS;
  }

  static HiddenSemiMarkov*
  hidden_semi_markov_stochastic_estimation_model(
      const MarkovianSequences &input, char model_type, int nb_state,
      bool left_right, double occupancy_mean, bool common_dispersion, int min_nb_state_sequence,
      int max_nb_state_sequence, double parameter, int estimator, bool counting_flag,
      bool state_sequence, int nb_iter)
  {
    HEADER_OS(HiddenSemiMarkov);
    ret = input.hidden_semi_markov_stochastic_estimation(error, os, model_type,
        nb_state, left_right, occupancy_mean, common_dispersion, min_nb_state_sequence,
        max_nb_state_sequence, parameter, estimator, counting_flag, state_sequence,
        nb_iter);
    FOOTER_OS;

  }

  static SemiMarkov*
  semi_markov_estimation(const MarkovianSequences &input, char model_type,
      int estimator, bool counting_flag, int nb_iter, int mean_computation_method)
  {
    HEADER_OS(SemiMarkov);
    ret = input.semi_markov_estimation(error, os, model_type, estimator,
        counting_flag, nb_iter, mean_computation_method);
    FOOTER_OS;
  }

  static bool
  mtg_write(const MarkovianSequences &input, const char* path,
      boost::python::list input_values)
  {
    StatError error;
    bool ret;

    CREATE_ARRAY(input_values, int, data);
    ret = input.mtg_write(error, path, data.get());

    return ret;
  }

  static NonhomogeneousMarkov*
  nonhomogeneous_markov_estimation(const MarkovianSequences &input,
      boost::python::list list_ident, bool counting_flag)
  {
    HEADER(NonhomogeneousMarkov);
    CREATE_ARRAY(list_ident, int, data);
    ret = input.nonhomogeneous_markov_estimation(error, data.get(), counting_flag);
    FOOTER;
  }

  static bool
  lumpability_test(const MarkovianSequences &input,
      boost::python::list input_values, int order)
  {
    StatError error;
    bool ret;
    std::stringstream os;

    CREATE_ARRAY(input_values, int, data);
    ret = input.lumpability_test(error, os, data.get(), order);
    FOOTER_OS;
  }

  static void
  plot_write(const MarkovianSequences &input, const std::string& prefix,
      const std::string& title)
  {
    StatError error;
    input.plot_write(error, prefix.c_str(), title.c_str());
  }

  static void
  file_ascii_sequences_write(const MarkovianSequences &input, const std::string& filename,
			     const std::string& format, bool exhaustive)
  {
    StatError error;
    if (format.length() > 0)
      input.ascii_data_write(error, filename.c_str(), format[0], exhaustive) ;
  }

  static SequenceCharacteristics *
  get_characteristics(const MarkovianSequences &input, const int variable)
  {
    // todo : shall we return something ? 
    HEADER(SequenceCharacteristics);
    ret = input.get_characteristics(variable);
    FOOTER;
  }

  static MultiPlotSet*
  get_plotable(const MarkovianSequences &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }

  static MarkovianSequences* 
  initial_run_computation(const MarkovianSequences &input)
  {
    HEADER(MarkovianSequences);
    ret = input.initial_run_computation(error);
    FOOTER;
  }


};

// Boost declaration

void class_markovian_sequences() {




  class_<MarkovianSequences, bases<Sequences> > ("_MarkovianSequences", "MarkovianSequences")

     //todo those constructors if needed
    .def("__init__", make_constructor(WRAP::markovian_sequences1))
    .def("__init__", make_constructor(WRAP::markovian_sequences2))
    .def("__init__", make_constructor(WRAP::markovian_sequences3))

    .def(init <Sequences>())

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
    .def("file_ascii_sequences_write", WRAP::file_ascii_sequences_write, "file name : string, format : character, exhaustive : boolean")
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
    //DEF_RETURN_VALUE("mtg_write", WRAP::mtg_write, args(""), "mtg_write")

    .def("mtg_write", WRAP::mtg_write, args(""), "")
    .def("plot_write", WRAP::plot_write, args("prefix", "title"), "Write GNUPLOT files")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable,"get_plotable")
    DEF_RETURN_VALUE("get_characteristics", WRAP::get_characteristics, args("variable"), "get_plotable")
    DEF_RETURN_VALUE("initial_run_computation", WRAP::initial_run_computation, args(""), "initial run computation")



    ;


  /*
   MarkovianSequences* remove_variable_1() const;
   MarkovianSequences* initial_run_computation(StatError &error) const;

   std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,   bool exhaustive = false) const;
   bool ascii_data_write(StatError &error , const char *path , char format = 'c' , bool exhaustive = false) const;
   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(StatError &error , const char *path ,   bool exhaustive = false) const;
   bool spreadsheet_write(StatError &error , const char *path) const;
   bool plot_write(StatError &error , const char *prefix ,   const char *title = 0) const;

   bool mtg_write(StatError &error , const char *path , int *itype) const;

   double iid_information_computation() const;

   void sojourn_time_histogram_computation(int variable);
   void create_observation_histogram(int nb_state);
   void observation_histogram_computation();
   void build_observation_histogram();
   void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,  bool initial_run_flag = false);

   // acces membres de la classe

   Curves* get_self_transition(int state) const { return self_transition[state]; }
   FrequencyDistribution*** get_observation() const { return observation; }
   FrequencyDistribution** get_observation(int variable) const { return observation[variable]; }
   FrequencyDistribution* get_observation(int variable , int state) const   { return observation[variable][state]; }
   */

}

#undef WRAP


void class_self_transition() {

  class_<SelfTransition, bases<Curves> > ("_SelfTransition", "SelfTransition", no_init)
    .def(init <int>())
    DEF_RETURN_VALUE_NO_ARGS("monomolecular_regression", &SelfTransition::monomolecular_regression, "returns monomolecular regression")
    DEF_RETURN_VALUE_NO_ARGS("logistic_regression", &SelfTransition::logistic_regression,"returns logistic regression")
    ;
    //Done
}

