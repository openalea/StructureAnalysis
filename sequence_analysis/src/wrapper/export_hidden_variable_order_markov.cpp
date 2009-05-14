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
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "stat_tool/distance_matrix.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
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
//using namespace stat_tool;

#define WRAP HiddenVariableOrderMarkovWrap
class HiddenVariableOrderMarkovWrap {

public:

  static boost::shared_ptr<Hidden_variable_order_markov>
  hidden_variable_order_markov_from_file(char* filename, int length,
		  double cumul_threshold)
  {
    Format_error error;

	Hidden_variable_order_markov *hvom = NULL;

	hvom = hidden_variable_order_markov_ascii_read(error, filename,
			length,  cumul_threshold);

	return boost::shared_ptr<Hidden_variable_order_markov>(hvom);
  }

  static Hidden_variable_order_markov*
  thresholding(const Hidden_variable_order_markov& input, double min_probability = MIN_PROBABILITY)
  {
	  //todo check this function
	Hidden_variable_order_markov *ret = NULL;
    ret =  input.thresholding(min_probability);
    return ret;
  }

  static Variable_order_markov_data*
  simulation_markovian_sequences(const Hidden_variable_order_markov &input, int nb_sequence,
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
  state_sequence_computation(const Hidden_variable_order_markov& input,
  		  const Markovian_sequences &seq , bool characteristic_flag = true)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, state_sequence_computation,
  		Variable_order_markov_data,  seq, characteristic_flag);
  }

  static Variable_order_markov_data*
  simulation_histogram(const Hidden_variable_order_markov &input,
		  const Histogram &hlength,  bool counting_flag, bool divergence_flag)
  {
	Format_error error;
	Variable_order_markov_data* ret;

	ret = input.simulation(error, hlength, counting_flag, divergence_flag);

	return ret;
  }

  static Variable_order_markov_data*
  simulation_nb_sequences(const Hidden_variable_order_markov &input,
		  int nb_sequence, int length, bool counting_flag)
  {
	Format_error error;
	Variable_order_markov_data* ret;

	ret = input.simulation(error,nb_sequence, length, counting_flag);

	return ret;
  }


  //Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,st Hidden_variable_order_markov **ihmarkov , Histogram **hlength , const char *path = 0) const;

  //Distance_matrix* divergence_computation(Format_error &error , std::ostream &os ,int nb_model ,const Hidden_variable_order_markov **hmarkov , int nb_sequence ,   int length , const char *path = 0) const;



  static Distance_matrix*
  divergence_computation(const Hidden_variable_order_markov &input,
      boost::python::list &input_markov, boost::python::list &input_sequence,
      int nb_seq, char *filename)
  {
    // there is as much Variable_order_markov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    Distance_matrix *ret = NULL;
    Format_error error;
    std::stringstream os;
    int nb_markov = len(input_markov);
    int nb_sequence = len(input_sequence);

    sequence_analysis::wrap_util::auto_ptr_array<
        const Hidden_variable_order_markov *> markov(
        new const Hidden_variable_order_markov*[nb_markov]);
    for (int i = 0; i < nb_markov; i++)
      markov[i] = boost::python::extract<Hidden_variable_order_markov*>(
          input_markov[i]);

    sequence_analysis::wrap_util::auto_ptr_array<const Markovian_sequences *>
        sequence(new const Markovian_sequences*[nb_sequence]);
    for (int i = 0; i < nb_sequence; i++)
      sequence[i] = boost::python::extract<Markovian_sequences*>(
          input_sequence[i]);

    ret = input.divergence_computation(error, os, nb_markov + 1, markov.get(),
        nb_seq, sequence.get(), filename);
    cerr << os.str() << endl;
    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);

    return ret;
  }
};



void class_hidden_variable_order_markov() {

        class_<Hidden_variable_order_markov, bases<Variable_order_markov > >
        ("_Hidden_variable_order_markov", "Hidden_variable_order_markov")
        .def("__init__", make_constructor(WRAP::hidden_variable_order_markov_from_file))


        .def(self_ns::str(self)) //__str__

		DEF_RETURN_VALUE("thresholding", WRAP::thresholding, args("probability"), "todo")
		DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
		DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
		DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
		DEF_RETURN_VALUE("state_sequence_computation", WRAP::state_sequence_computation, args(""),"")
		DEF_RETURN_VALUE("divergence_computation", WRAP::divergence_computation, args("input", "input_markov", "input_sequence", "filename"), "todo")

        ;
      /*
        Hidden_variable_order_markov() {}
          Hidden_variable_order_markov(const Variable_order_markov *pmarkov , int inb_output_process ,
                                       Nonparametric_process **nonparametric_observation ,
                                       Parametric_process **parametric_observation , int length)

          Hidden_variable_order_markov(const Hidden_variable_order_markov &hmarkov , bool data_flag = true)

          std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
          bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
          bool spreadsheet_write(Format_error &error , const char *path) const;

		 double likelihood_computation(const Markovian_sequences &seq ,double *posterior_probability = 0 ,  int index = I_DEFAULT) const;

          bool state_profile_write(Format_error &error , std::ostream &os , const Markovian_sequences &iseq ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,            int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(Format_error &error , const char *path , const Markovian_sequences &iseq ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_ascii_write(Format_error &error , std::ostream &os , int identifier ,int state_sequence = GENERALIZED_VITERBI ,int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(Format_error &error , const char *path ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,    int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_plot_write(Format_error &error , const char *prefix ,  const Markovian_sequences &iseq ,    int identifier , const char *title = 0) const;
          bool state_profile_plot_write(Format_error &error , const char *prefix ,   int identifier , const char *title = 0) const;




*/
}


#undef WRAP
