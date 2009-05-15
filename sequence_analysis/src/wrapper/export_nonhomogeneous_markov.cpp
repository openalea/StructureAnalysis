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
#include "stat_tool/regression.h"
#include "sequence_analysis/sequences.h"
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




#define WRAP NonHomogeneousMarkovWrap
class NonHomogeneousMarkovWrap
{

public:


  static boost::shared_ptr<Nonhomogeneous_markov>
    nonhomogeneous_markov_from_file(char* filename, int length)
    {
      //olf_format should be true
      Format_error error;
      Nonhomogeneous_markov *nonhomo = NULL;
      nonhomo = nonhomogeneous_markov_ascii_read(error,
          filename, length);
      if(!nonhomo)
      {
        sequence_analysis::wrap_util::throw_error(error);
      }
      return boost::shared_ptr<Nonhomogeneous_markov>(nonhomo);
    }

  static Nonhomogeneous_markov_data*
  simulation_markovian_sequences(const Nonhomogeneous_markov &input,
      int nb_sequence, const Markovian_sequences input_seq, bool counting_flag)
  {
    Format_error error;
    Nonhomogeneous_markov_data* ret = NULL;

    ret = input.simulation(error, nb_sequence, input_seq, counting_flag);

    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);

    return ret;
  }

  static Nonhomogeneous_markov_data*
  simulation_histogram(const Nonhomogeneous_markov &input,
      const Histogram &hlength, bool counting_flag)
  {
    Format_error error;
    Nonhomogeneous_markov_data* ret;

    ret = input.simulation(error, hlength, counting_flag);

    return ret;
  }

  static Nonhomogeneous_markov_data*
  simulation_nb_sequences(const Nonhomogeneous_markov &input, int nb_sequence,
      int length, bool counting_flag)
  {
    Format_error error;
    Nonhomogeneous_markov_data* ret;

    ret = input.simulation(error, nb_sequence, length, counting_flag);

    return ret;
  }

  static Parametric_model*
    extract(const Nonhomogeneous_markov& seq, int type, int state)
    {
      Format_error error;
      Parametric_model* ret;
      ret = seq.extract(error, type, state);
      if (!ret)
        sequence_analysis::wrap_util::throw_error(error);
      return ret;
    }

};

// Boost declaration

void class_nonhomogeneous_markov() {

  class_<Nonhomogeneous_markov, bases<STAT_interface> > ("_Nonhomogeneous_markov", "Nonhomogeneous_markov")
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_nb_state_and_ident_list))
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_chain_and_self_transition))
    .def("__init__", make_constructor(WRAP::nonhomogeneous_markov_from_file))

    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")
    .def("get_homogeneity", &Nonhomogeneous_markov::get_homogeneity, python::args("index"),"return homogeneity")
    .def("get_self_transition", &Nonhomogeneous_markov::get_self_transition, return_value_policy<manage_new_object> (),	python::args("index"),"return Function* of self transition")

    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")

	//.def("get_process", &Nonhomogeneous_markov::get_process,
	//	"return non parametric sequence process")

/*
    Nonhomogeneous_markov();
    Nonhomogeneous_markov(int inb_state , int *ident);
    Nonhomogeneous_markov(const Chain *pchain , const Function **pself_transition , int length);
    Nonhomogeneous_markov(const Nonhomogeneous_markov &markov , bool data_flag = true , bool haracteristic_flag = true)    :Chain(markov) { copy(markov , data_flag , characteristic_flag); }


    void characteristic_computation(int length , bool counting_flag);
    void characteristic_computation(const Nonhomogeneous_markov_data &seq , bool counting_flag ,            bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq , int index = I_DEFAULT) const;

    Nonhomogeneous_markov_data* get_markov_data() const { return markov_data; }
    Nonparametric_sequence_process* get_process() const { return process; }
   */



;
}
#undef WRAP



#define WRAP NonHomogeneousMarkovDataWrap
class NonHomogeneousMarkovDataWrap
{

public:
  static Distribution_data*
  extract(const Nonhomogeneous_markov_data& seq, int type, int state)
  {
    Format_error error;
    Distribution_data* ret;
    ret = seq.extract(error, type, state);
    if (!ret)
      sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }

  static Nonhomogeneous_markov_data*
  remove_index_parameter(const Nonhomogeneous_markov_data& seq)
  {
    SIMPLE_METHOD_TEMPLATE_0(seq, remove_index_parameter, Nonhomogeneous_markov_data);
  }

};

void class_nonhomogeneous_markov_data()
{
  class_<Nonhomogeneous_markov_data, bases<Markovian_sequences > >
  ("_Nonhomogeneous_markov_data", "Nonhomogeneous_markov_data")

    .def(init<Markovian_sequences & >())

    .add_property("likelihood", &Nonhomogeneous_markov_data::get_likelihood, "returns likelihood")
    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")

  ;
  /*
    Nonhomogeneous_markov_data();
    Nonhomogeneous_markov_data(const Histogram &ihlength);
    Nonhomogeneous_markov_data(const Nonhomogeneous_markov_data &seq , bool model_flag = true , char transform = 'c') :Markovian_sequences(seq , transform) { copy(seq , model_flag); }

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,const char *title = 0) const;

    void build_transition_count();

    Nonhomogeneous_markov* get_markov() const { return markov; }
    Chain_data* get_chain_data() const { return chain_data; }

    */

}

#undef WRAP
