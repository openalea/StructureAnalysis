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
 *        $Id: export_tops.cpp 6169 2009-04-01 16:42:59Z cokelaer $
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
#include "sequence_analysis/variable_order_markov.h"
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

class MarkovianSequencesWrap {

public:



  static Markovian_sequences*
  markovian_sequences1(int nb_sequence, boost::python::list& input_identifiers,
      boost::python::list& input_ilength, int input_index_parameter_type,
      int input_nb_variable, boost::python::list& input_itype, bool init_flag)
  {
    Markovian_sequences *seq = NULL;
    // todo
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

  // Merge
  static Markovian_sequences*
  merge(const Markovian_sequences& input, const boost::python::list& input_list)
  {
    int size = len(input_list);
    sequence_analysis::wrap_util::auto_ptr_array<const Markovian_sequences *>
        sequens(new const Markovian_sequences*[size]);
    for (int i = 0; i < size; i++)
      sequens[i] = boost::python::extract<Markovian_sequences *>(input_list[i]);

    SIMPLE_METHOD_TEMPLATE_1(input, merge, Markovian_sequences, size, sequens.get());
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
		  int variable, boost::python::list& limit, bool add_flag = false)
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
           ldouble[i] = boost::python::extract<int> (limit[i]);
         else
           lint[i] = boost::python::extract<double> (limit[i]);
       }

     // Call correct function
     if (is_float)
       {
         ret = seq.cluster(error, variable, nb_limit, ldouble);
         delete[] ldouble;
       }
     else
       {
         ret = seq.cluster(error, variable, nb_limit, lint, add_flag);
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

    CREATE_ARRAY(input_list, int)
    //int size = len(input_list);
    //sequence_analysis::wrap_util::auto_ptr_array<int>
    //  l(new int[size]);
    //for (int i = 0; i < size; i++)
    //  l[i] = boost::python::extract<int>(input_list[i]);

    SIMPLE_METHOD_TEMPLATE_1(input, transcode, Markovian_sequences,
    		variable, data.get(), add_flag);
  }


  static Markovian_sequences*
  select_variable(const Markovian_sequences& input, int inb_variable,
      boost::python::list& input_list , bool keep = true)
  {
    int size = len(input_list);
    sequence_analysis::wrap_util::auto_ptr_array<int> l(new int[size]);
    for (int i = 0; i < size; i++)
      l[i] = boost::python::extract<int>(input_list[i]);

    SIMPLE_METHOD_TEMPLATE_1(input, select_variable, Markovian_sequences, inb_variable, l.get(), keep);
  }

  static Markovian_sequences*
  split(const Markovian_sequences& input, int index)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, split, Markovian_sequences, index);
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
	  //sequence_length = I_DEFAULT , int run_length = I_DEFAULT)
     SIMPLE_METHOD_TEMPLATE_1(input, add_absorbing_run,
    		 Markovian_sequences, sequence_length, run_length);
  }


  static std::string
  word_count(const Markovian_sequences& input, int variable, int word_length,
		  int begin_state, int end_state, int min_frequency)
  {
	  //default values for begin_state=I_DEFAULT
	  //default values for end_state=I_DEFAULT
	  //default values for min_frequency=1
	Format_error error;
	std::stringstream s;


	bool res = true;
	res = input.word_count(error, s ,variable, word_length,
			begin_state, end_state,  min_frequency);
	if (!res)
		sequence_analysis::wrap_util::throw_error(error);

	return s.str();
  }


  static Variable_order_markov*
  variable_order_markov_estimation1(const Markovian_sequences& input,
		  char model_type, int min_order, int max_order, int algorithm,
		  double threshold, int estimator, bool global_initial_transition,
		  bool global_sample, bool counting_flag)
  {
	  Format_error error;
	  Variable_order_markov *vom;
	  std::stringstream os;

	  vom = input.variable_order_markov_estimation(error, os, model_type, min_order,
			  max_order, algorithm, threshold, estimator, global_initial_transition,
			  global_sample, counting_flag);

	  if (!vom)
	  		sequence_analysis::wrap_util::throw_error(error);
	  return vom;
  }

  static Variable_order_markov*
  variable_order_markov_estimation2(const Markovian_sequences& input,
		  char type, int max_order, bool global_initial_transition,
		  bool counting_flag)
  {
	  Format_error error;
	  Variable_order_markov *vom;
	  std::stringstream os;

	  vom = input.variable_order_markov_estimation(error, type,
			  max_order, global_initial_transition, counting_flag);

	  if (!vom)
	  		sequence_analysis::wrap_util::throw_error(error);
	  return vom;
  }

  static Variable_order_markov*
  variable_order_markov_estimation3(const Markovian_sequences& input,
		  const Variable_order_markov &markov,
		  bool global_initial_transition,
		  bool counting_flag)
  {
	  Format_error error;
	  Variable_order_markov *vom;
	  std::stringstream os;

	  vom = input.variable_order_markov_estimation(error, markov,
			  global_initial_transition, counting_flag);

	  if (!vom)
	  		sequence_analysis::wrap_util::throw_error(error);
	  return vom;
  }

  static Variable_order_markov*
  lumpability_estimation(const Markovian_sequences& input,
		  boost::python::list& input_symbol,
		  int penalty_type,
		  int order,
		  bool counting_flag)
  {
    Format_error error;
    Variable_order_markov *vom;
    std::stringstream os;

    int nb_value = len(input_symbol);
    int *symbol;
    symbol = new int[nb_value];

    for (int i = 0; i < nb_value; i++)
    {
      symbol[i] = boost::python::extract<int> (input_symbol[i]);
    }

    vom = input.lumpability_estimation(error, os, symbol,
		  penalty_type, order, counting_flag);

    if (!vom)
  		sequence_analysis::wrap_util::throw_error(error);

    delete[] symbol;
    return vom;
  }

};

// Boost declaration

void class_markovian_sequences() {

  enum_<sequence_analysis::wrap_util::UniqueInt<11, 102> >("MarkovianSequenceType")
      .value("OBSERVATION",OBSERVATION)
      .value("FIRST_OCCURRENCE",FIRST_OCCURRENCE)
      .value("RECURRENCE_TIME",RECURRENCE_TIME)
      .value("SOJOURN_TIME",SOJOURN_TIME)
      .value("INITIAL_RUN",INITIAL_RUN)
      .value("FINAL_RUN",FINAL_RUN)
      .value("NB_RUN",NB_RUN)
      .value("NB_OCCURRENCE",NB_OCCURRENCE)
      .value("LENGTH",LENGTH)
      .value("SEQUENCE_CUMUL",SEQUENCE_CUMUL)
      .value("SEQUENCE_MEAN",SEQUENCE_MEAN)
      .export_values();

  class_<Markovian_sequences, bases<Sequences> > ("_Markovian_sequences", "Markovian_sequences")

    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences1))
    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences2))
    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences3))
    .def(init <Sequences>())
    .def(init <Markovian_sequences, optional<char, int> >())

    .def(self_ns::str(self)) //__str__

    DEF_RETURN_VALUE("extract", MarkovianSequencesWrap::extract,args("type", "variable","value"), "Extract distribution data")
    DEF_RETURN_VALUE("split", &MarkovianSequencesWrap::split,args("step"),  "Split")
    DEF_RETURN_VALUE("merge", &MarkovianSequencesWrap::merge,args("sequences"),  "Merge")
    DEF_RETURN_VALUE("cluster_step", &MarkovianSequencesWrap::cluster_step,args("variable", "step", "mode"), "Cluster")
    DEF_RETURN_VALUE("cluster_limit", &MarkovianSequencesWrap::cluster_limit,args("variable", "limit", "bool_add_flag"), "Cluster")
    DEF_RETURN_VALUE("transcode", &MarkovianSequencesWrap::transcode,args("variable", "symbol", "add_flag"), "Transcode")
    DEF_RETURN_VALUE("select_variable", &MarkovianSequencesWrap::select_variable,args("nb_variable", "list variables", "keep"), "select variable")
    DEF_RETURN_VALUE("add_absorbing_run", &MarkovianSequencesWrap::add_absorbing_run,args("sequence_length", "run_length"), "todo")
    .def("word_count", &MarkovianSequencesWrap::word_count, args("variable", "word_length","begin_state", "end_state","min_frequency"), "todo" )
    DEF_RETURN_VALUE_NO_ARGS("remove_index_parameter", &MarkovianSequencesWrap::remove_index_parameter, "Remove index parameter")




    DEF_RETURN_VALUE("variable_order_markov_estimation1", &MarkovianSequencesWrap::variable_order_markov_estimation1, args("model_type", "min_order", "max_order", "algorithm", "threshold", "estimator","global_initial_transition","global_sample", "counting_flag"), "todo")
    DEF_RETURN_VALUE("variable_order_markov_estimation2", &WRAP::variable_order_markov_estimation2, args("type","max_order","global_initial_transition","counting_flag"), "todo")
    DEF_RETURN_VALUE("variable_order_markov_estimation3",&WRAP::variable_order_markov_estimation3,args("markov","global_initial_transition","counting_flag"), "todo")
    DEF_RETURN_VALUE("lumpability_estimation", &WRAP::lumpability_estimation, args("input_symbol", "penalty_type","order","counting_flag"), "todo")//

    ;
	/*


   Markovian_sequences* consecutive_values(Format_error &error , std::ostream &os , int ivariable , bool add_flag = false) const;
   Markovian_sequences* cluster(Format_error &error , int ivariable , int nb_class , int *ilimit , bool add_flag = false) const;
   Markovian_sequences* cluster(Format_error &error , int variable , int nb_class , double *ilimit) const;


   Markovian_sequences* merge_variable(Format_error &error , int nb_sample , const Markovian_sequences **iseq , int ref_sample = I_DEFAULT) const;
   Markovian_sequences* remove_variable_1() const;

   Markovian_sequences* initial_run_computation(Format_error &error) const;


   std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,   bool exhaustive = false) const;
   bool ascii_data_write(Format_error &error , const char *path , char format = 'c' , bool exhaustive = false) const;
   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(Format_error &error , const char *path ,   bool exhaustive = false) const;
   bool spreadsheet_write(Format_error &error , const char *path) const;
   bool plot_write(Format_error &error , const char *prefix ,   const char *title = 0) const;

   bool transition_count(Format_error &error , std::ostream &os , int max_order , bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,  const char *path = 0) const;
   bool mtg_write(Format_error &error , const char *path , int *itype) const;

   double iid_information_computation() const;

   void self_transition_computation();
   void self_transition_computation(bool *homogeneity);
   void sojourn_time_histogram_computation(int variable);
   void create_observation_histogram(int nb_state);
   void observation_histogram_computation();
   void build_observation_histogram();
   void build_characteristic(int variable = I_DEFAULT , bool sojourn_time_flag = true ,  bool initial_run_flag = false);

   Nonhomogeneous_markov* nonhomogeneous_markov_estimation(Format_error &error , int *ident ,  bool counting_flag = true) const;







   Semi_markov* semi_markov_estimation(Format_error &error , std::ostream &os , char model_type ,
                                       int estimator = COMPLETE_LIKELIHOOD , bool counting_flag = true ,
                                       int nb_iter = I_DEFAULT , int mean_computation = COMPUTED) const;

   Hidden_variable_order_markov* hidden_variable_order_markov_estimation(Format_error &error , std::ostream &os ,
                                                                         const Hidden_variable_order_markov &ihmarkov ,
                                                                         bool global_initial_transition = true ,
                                                                         bool counting_flag = true ,
                                                                         bool state_sequence = true ,
                                                                         int nb_iter = I_DEFAULT) const;
   Hidden_variable_order_markov* hidden_variable_order_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                                    const Hidden_variable_order_markov &ihmarkov ,
                                                                                    bool global_initial_transition = true ,
                                                                                    int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                                    int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                                    double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                                    bool counting_flag = true ,
                                                                                    bool state_sequence = true ,
                                                                                    int nb_iter = I_DEFAULT) const;

   Hidden_semi_markov* hidden_semi_markov_estimation(Format_error &error , std::ostream &os ,
                                                     const Hidden_semi_markov &ihsmarkov ,
                                                     int estimator = COMPLETE_LIKELIHOOD ,
                                                     bool counting_flag = true ,
                                                     bool state_sequence = true ,
                                                     int nb_iter = I_DEFAULT ,
                                                     int mean_computation = COMPUTED) const;
   Hidden_semi_markov* hidden_semi_markov_estimation(Format_error &error , std::ostream &os ,
                                                     char model_type , int nb_state , bool left_right ,
                                                     int estimator = COMPLETE_LIKELIHOOD ,
                                                     bool counting_flag = true ,
                                                     bool state_sequence = true ,
                                                     double occupancy_mean = D_DEFAULT ,
                                                     int nb_iter = I_DEFAULT ,
                                                     int mean_computation = COMPUTED) const;
   Hidden_semi_markov* hidden_semi_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                const Hidden_semi_markov &ihsmarkov ,
                                                                int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                int estimator = COMPLETE_LIKELIHOOD ,
                                                                bool counting_flag = true ,
                                                                bool state_sequence = true ,
                                                                int nb_iter = I_DEFAULT) const;
   Hidden_semi_markov* hidden_semi_markov_stochastic_estimation(Format_error &error , std::ostream &os ,
                                                                char model_type , int nb_state , bool left_right ,
                                                                int min_nb_state_sequence = MIN_NB_STATE_SEQUENCE ,
                                                                int max_nb_state_sequence = MAX_NB_STATE_SEQUENCE ,
                                                                double parameter = NB_STATE_SEQUENCE_PARAMETER ,
                                                                int estimator = COMPLETE_LIKELIHOOD ,
                                                                bool counting_flag = true ,
                                                                bool state_sequence = true ,
                                                                double occupancy_mean = D_DEFAULT ,
                                                                int nb_iter = I_DEFAULT) const;

   bool lumpability_test(Format_error &error , std::ostream &os , int *symbol , int order = 1) const;
   bool comparison(Format_error &error , std::ostream &os , int nb_model ,   const Variable_order_markov **imarkov , const char *path = 0) const;
   bool comparison(Format_error &error , std::ostream &os , int nb_model , const Semi_markov **ismarkov , const char *path = 0) const
   bool comparison(Format_error &error , std::ostream &os , int nb_model ,  const Hidden_variable_order_markov **ihmarkov , int algorithm = FORWARD , const char *path = 0) const;
   bool comparison(Format_error &error , std::ostream &os , int nb_model ,  const Hidden_semi_markov **ihsmarkov ,int algorithm = FORWARD , const char *path = 0) const;

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

    DEF_RETURN_VALUE_NO_ARGS("monomolecular_regression", &Self_transition::monomolecular_regression, "return function monomolecular regression")
    DEF_RETURN_VALUE_NO_ARGS("logistic_regression", &Self_transition::logistic_regression,"return function logistic regression")
    ;
    //Done
}
