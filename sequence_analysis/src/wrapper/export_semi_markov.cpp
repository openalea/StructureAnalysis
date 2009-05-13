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

    /*    if(!top_parameters)
     {
     sequence_analysis::wrap_util::throw_error(error);
     }
     */
    return boost::shared_ptr<Semi_markov>(semi_markov);
  }

  static Semi_markov_data*
  extact_data(const Semi_markov& input)
  {
    Format_error err;
    Semi_markov_data *ret;
    ret = input.extract_data(err);
    return ret;
    //SIMPLE_METHOD_TEMPLATE_1(input, extract_data, Semi_markov_data);
  }

  static Parametric_model*
  extract_histogram(const Semi_markov& input, int state, int histogram_type = FINAL_RUN)
  {
    // does not work issue with sojourn...
    SIMPLE_METHOD_TEMPLATE_1(input, extract, Parametric_model, state, histogram_type);
  }
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(extract_histogram_overloads, SemiMarkovWrap::extract_histogram, 1, 2);


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




  /*
  static Time_events*
    nb_event_select(const Time_events& input, int min, int max){
      SIMPLE_METHOD_TEMPLATE_1(input, nb_event_select, Time_events, min, max);
    }
*/
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
    .def("extract_histogram", (Parametric_model *(*)(const Semi_markov&, int, int))  SemiMarkovWrap::extract_histogram, return_value_policy< manage_new_object >(), SemiMarkovWrap::extract_histogram_overloads())
    .def("extract_histogram", (Parametric_model *(*)(const Semi_markov& , int))  SemiMarkovWrap::extract_histogram, return_value_policy< manage_new_object >(), SemiMarkovWrap::extract_histogram_overloads())

    .def("get_forward", (Forward *(*)(const Semi_markov&, int)) SemiMarkovWrap::get_forward, return_value_policy< manage_new_object >(), SemiMarkovWrap::get_forward_overloads())

    DEF_RETURN_VALUE_NO_ARGS("get_semi_markov_data", &Semi_markov::get_semi_markov_data, "returns semi_markov_data")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", &Semi_markov::extract_data, "returns semi_markov_data")

    .def("file_ascii_write", SemiMarkovWrap::file_ascii_write,"Save vector summary into a file")

    ;


  //todo file_ascii_write not accessible ?
  /*
   *
   Semi_markov();
   Semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value);
   Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy , const Nonparametric_process *pobservation , int length ,  bool counting_flag);
   Semi_markov(const Semi_markov &smarkov , bool data_flag = true ,  int param = I_DEFAULT)  :Chain(smarkov) { copy(smarkov , data_flag , param); }

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;


   Semi_markov* thresholding(double min_probability = MIN_PROBABILITY) const;
   void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
   void characteristic_computation(const Semi_markov_data &seq , bool counting_flag, int variable = I_DEFAULT , bool length_flag = true);

   double likelihood_computation(const Markovian_sequences &seq , int index) const;
   double likelihood_computation(const Semi_markov_data &seq) const;

   Semi_markov_data* simulation(Format_error &error , const Histogram &hlength ,   bool counting_flag = true , bool divergence_flag = false) const;
   Semi_markov_data* simulation(Format_error &error , int nb_sequence , int length , bool counting_flag = true) const;
   Semi_markov_data* simulation(Format_error &error , int nb_sequence , const Markovian_sequences &iseq , bool counting_flag = true) const;

   Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,   const Semi_markov **ismarkov , Histogram **hlength ,   const char *path = 0) const;
   Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,   const Semi_markov **smarkov , int nb_sequence ,   int length , const char *path = 0) const;
   Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,   const Semi_markov **smarkov , int nb_sequence ,   const Markovian_sequences **seq , const char *path = 0) const;

   Forward** get_forward() const { return forward; }

   Nonparametric_sequence_process* get_nonparametric_process(int variable)   const { return nonparametric_process[variable]; }
   Parametric_process** get_parametric_process() const { return parametric_process; }


   */

  ;
}


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

