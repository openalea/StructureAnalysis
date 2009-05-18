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
#include "stat_tool/stat_label.h"
#include "stat_tool/distance_matrix.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
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


#define WRAP VariableOrderMarkovWrap
class VariableOrderMarkovWrap {

public:

  static boost::shared_ptr<Variable_order_markov>
  variable_order_markov_from_file(char* filename)
  {
    Format_error error;
    Variable_order_markov *vom = NULL;
    vom = variable_order_markov_ascii_read(error, filename);
    if(!vom)
    {
      sequence_analysis::wrap_util::throw_error(error);
    }
    return boost::shared_ptr<Variable_order_markov>(vom);
  }


  static void
   file_ascii_write(const Variable_order_markov& d, const char* path, bool exhaustive)
   {
     bool result = true;
     Format_error error;

     result = d.ascii_write(error, path, exhaustive);
     if (!result)
       sequence_analysis::wrap_util::throw_error(error);
   }

  static Parametric_model*
  extract(const Variable_order_markov& input, int type , int variable , int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, Parametric_model, type, variable, value);
  }

  static Variable_order_markov_data*
  extract_data(const Variable_order_markov& input)
  {
    SIMPLE_METHOD_TEMPLATE_0(input, extract_data, Variable_order_markov_data);
  }


  static Correlation*
  output_autocorrelation_computation(const Variable_order_markov& input,
		  int variable, int output, int max_lag)
  {
	 SIMPLE_METHOD_TEMPLATE_1(input, output_autocorrelation_computation,
			 Correlation, variable, output, max_lag);
  }

  static Correlation*
  state_autocorrelation_computation(const Variable_order_markov& input,
		  int state, int max_lag)
  {
	 SIMPLE_METHOD_TEMPLATE_1(input, state_autocorrelation_computation,
			 Correlation, state, max_lag);
  }


  static double
  likelihood_computation(const Variable_order_markov& input, const Markovian_sequences& ms, int index)
  {
	return input.likelihood_computation(ms, index);
  }

  static Variable_order_markov*
  thresholding(const Variable_order_markov& input, double min_probability = MIN_PROBABILITY)
  {
	return input.thresholding(min_probability);
  }

  static Distance_matrix*
  divergence_computation(const Variable_order_markov &input,
		  boost::python::list &input_markov,
		  boost::python::list &input_sequence,
		  int nb_seq,
		  char *filename)
  {
	// there is as much Variable_order_markov elmts as Markovian elts
	// 1 for input and N-1 for input_markov and N input_sequence

	Distance_matrix *ret = NULL;
	Format_error error;
    std::stringstream os;
    int nb_markov = len(input_markov);
    int nb_sequence = len(input_sequence);

    sequence_analysis::wrap_util::auto_ptr_array<const Variable_order_markov *>
		markov(new const Variable_order_markov*[nb_markov]);
    for (int i = 0; i < nb_markov; i++)
        markov[i] = boost::python::extract<Variable_order_markov*> (input_markov[i]);

    sequence_analysis::wrap_util::auto_ptr_array<const Markovian_sequences *>
   		sequence(new const Markovian_sequences*[nb_sequence]);
    for (int i = 0; i < nb_sequence; i++)
        sequence[i] = boost::python::extract<Markovian_sequences*> (input_sequence[i]);

	ret = input.divergence_computation(error, os, nb_markov+1, markov.get(),
			nb_seq, sequence.get(), filename);
    cerr << os.str() <<endl;
	if (!ret)
		sequence_analysis::wrap_util::throw_error(error);

	return ret;
  }


  /*Variable_order_markov_data* simulation(Format_error &error , const Histogram &hlength , bool counting_flag = true , bool divergence_flag = false) const;
  Variable_order_markov_data* simulation(Format_error &error , int nb_sequence , int length , bool counting_flag = true) const;
  Variable_order_markov_data* simulation(Format_error &error , int nb_sequence , const Markovian_sequences &iseq , bool counting_flag = true) const;
  */

  static Variable_order_markov_data*
    simulation_markovian_sequences(const Variable_order_markov &input, int nb_sequence,
  		  const Markovian_sequences input_seq, bool counting_flag)
    {
  	Format_error error;
  	Variable_order_markov_data* ret = NULL;

  	ret = input.simulation(error, nb_sequence,
  			input_seq, counting_flag);

      if (!ret)
          sequence_analysis::wrap_util::throw_error(error);

  	return ret;
    }
  static Variable_order_markov_data*
   simulation_histogram(const Variable_order_markov &input,
 		  const Histogram &hlength,  bool counting_flag, bool divergence_flag)
   {
 	Format_error error;
 	Variable_order_markov_data* ret;

 	ret = input.simulation(error, hlength, counting_flag, divergence_flag);

 	return ret;
   }

   static Variable_order_markov_data*
   simulation_nb_sequences(const Variable_order_markov &input,
 		  int nb_sequence, int length, bool counting_flag)
   {
 	Format_error error;
 	Variable_order_markov_data* ret;

 	ret = input.simulation(error,nb_sequence, length, counting_flag);

 	return ret;
   }


};

// Boost declaration

void class_variable_order_markov() {

  class_<Variable_order_markov, bases<STAT_interface > >
  ("_Variable_order_markov", "Variable_order_markov\n"
      "constructor(type(char in ['o','e']), nb_state(int), nb_row(int)\n"
      "constructor(type(char in ['o','e']), nb_state(int), nb_row(int), max_order(int)\n"
      "constructor(type(char in ['o','e']), nb_state(int), order(int), init_flag(bool), output_process=0, nb_value=0\n"

  )
    .def("__init__", make_constructor(VariableOrderMarkovWrap::variable_order_markov_from_file))



      // type = 'o' : ordinaire, or 'e' : en equilibre des probabilites de transition
    .def(init <char, int, int>())
    .def(init <char, int, int,int>())
    .def(init <char, int, int, bool, optional<int, int> >())

    .def(self_ns::str(self)) //__str__

    .def("file_ascii_write", WRAP::file_ascii_write,"Save vector summary into a file")
//    .def("ascii_data_write", WRAP::ascii_data_write,"Return a string with the object representation")
//    .def("file_ascii_data_write", WRAP::file_ascii_data_write,"Save vector data into a file")

    .add_property("nb_iterator", &Variable_order_markov::get_nb_iterator, "todo")
    .add_property("max_order", &Variable_order_markov::get_max_order, "todo")
    .add_property("nb_output_process", &Variable_order_markov::get_nb_output_process, "todo")

    .def("get_memory_type", &Variable_order_markov::get_memory_type, args("memory"),"todo")
    .def("get_order", &Variable_order_markov::get_order, args("memory"),"todo")
    .def("get_state", &Variable_order_markov::get_state, args("memory", "lag"),"todo")
    .def("get_parent", &Variable_order_markov::get_parent, args("memory"),"todo")
    .def("get_child", &Variable_order_markov::get_child, args("memory", "state"),"todo")
    .def("get_next", &Variable_order_markov::get_next, args("memory", "state"),"todo")
    .def("get_nb_memory", &Variable_order_markov::get_nb_memory, args("memory"),"todo")
    .def("get_previous", &Variable_order_markov::get_previous,args("memory", "state"), "todo")

    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data, "returns variable_order_markov_data")

    DEF_RETURN_VALUE("extract", WRAP::extract, args("type","variable","value"), "returns parametric model")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data, "returns variable_order_markov_data")
    DEF_RETURN_VALUE("output_autocorrelation_computation", WRAP::output_autocorrelation_computation, args("variable", "output", "max_lag"), "todo")
    DEF_RETURN_VALUE("state_autocorrelation_computation", WRAP::state_autocorrelation_computation, args("state", "max_lag"), "todo")

    .def("likelihood_computation", WRAP::likelihood_computation, args("seq", "index"),"todo")
    DEF_RETURN_VALUE("divergence_computation", WRAP::divergence_computation, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("thresholding", WRAP::thresholding, args("index"), "todo")

    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")


    ;


/*
 *
 *  friend Variable_order_markov* variable_order_markov_parsing(Format_error &error ,
                                                                ifstream &in_file ,
                                                                int &line , char type);

  Variable_order_markov(const Variable_order_markov &markov , int inb_output_process , int nb_value);
  Variable_order_markov(const Variable_order_markov *pmarkov ,  const Nonparametric_process *pobservation , int length);
  Variable_order_markov(const Variable_order_markov &markov , bool data_flag = true)  :Chain(markov) { copy(markov , data_flag); }

  void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
  void characteristic_computation(const Variable_order_markov_data &seq , bool counting_flag ,
  int variable = I_DEFAULT , bool length_flag = true);

  double likelihood_computation(const Variable_order_markov_data &seq) const;


  Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,const Variable_order_markov **imarkov , Histogram **hlength ,  const char *path = 0) const;
  Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,const Variable_order_markov **markov , int nb_sequence ,  int length , const char *path = 0) const;

  Variable_order_markov_data* get_markov_data() const { return markov_data; }
  Nonparametric_sequence_process* get_nonparametric_process(int variable) const{ return nonparametric_process[variable]; }
  Parametric_process** get_parametric_process() const { return parametric_process; }
  Parametric_process* get_parametric_process(int variable)const { return parametric_process[variable]; }
*/




}
#undef WRAP


void class_variable_order_markov_data() {

        class_<Variable_order_markov_data, bases<Markovian_sequences > >
        ("_Variable_order_markov_data", "Variable_order_markov_data")

;
        /*
Variable_order_markov_data();
   Variable_order_markov_data(const Histogram &ihlength , int inb_variable , bool init_flag = false);
   Variable_order_markov_data(const Markovian_sequences &seq);
   Variable_order_markov_data(const Markovian_sequences &seq , char transform ,
                              bool initial_run_flag);
   Variable_order_markov_data(const Variable_order_markov_data &seq , bool model_flag = true ,
                              char transform = 'c')
   :Markovian_sequences(seq , transform) { copy(seq , model_flag); }
   ~Variable_order_markov_data();
   Variable_order_markov_data& operator=(const Variable_order_markov_data &seq);

   Distribution_data* extract(Format_error &error , int type ,
                              int variable , int value) const;
   Variable_order_markov_data* remove_index_parameter(Format_error &error) const;

   Correlation* state_autocorrelation_computation(Format_error &error , int istate ,
                                                  int max_lag = MAX_LAG) const;
   Correlation* output_autocorrelation_computation(Format_error &error , int variable ,
                                                   int output , int max_lag = MAX_LAG) const;

   std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                  bool exhaustive = false) const;
   bool ascii_data_write(Format_error &error , const char *path ,
                         char format = 'c' , bool exhaustive = false) const;

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(Format_error &error , const char *path ,
                    bool exhaustive = false) const;
   bool spreadsheet_write(Format_error &error , const char *path) const;
   bool plot_write(Format_error &error , const char *prefix ,
                   const char *title = 0) const;

   void build_transition_count(const Variable_order_markov &markov ,
                               bool begin = true , bool non_terminal = false);
   void order0_estimation(Variable_order_markov &markov) const;

   // acces membres de la classe

   Variable_order_markov* get_markov() const { return markov; }
   Variable_order_chain_data* get_chain_data() const { return chain_data; }
   double get_likelihood() const { return likelihood; }
   double get_hidden_likelihood() const { return hidden_likelihood; }
   double get_posterior_probability(int index) const { return posterior_probability[index]; }
   */
}


