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
#include "stat_tool/stat_label.h"
#include "stat_tool/distance_matrix.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
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



#define WRAP HiddenSemiMarkovWrap
class WRAP {

public:

  static boost::shared_ptr<HiddenSemiMarkov>
  hidden_semi_markov_from_file(char *filename, int length, bool counting_flag,
      double cumul_threshold, bool old_format)
  {
    //olf_format should be true
    StatError error;
    HiddenSemiMarkov *hidden_semi_markov = NULL;

    hidden_semi_markov = hidden_semi_markov_ascii_read(error, filename, length,
        counting_flag, cumul_threshold, old_format);
    if (!hidden_semi_markov)
      {
        sequence_analysis::wrap_util::throw_error(error);
      }
    return boost::shared_ptr<HiddenSemiMarkov>(hidden_semi_markov);
  }

  static void
  file_ascii_write(const HiddenSemiMarkov &d, const char *path,
      bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }

  static SemiMarkovData*
  state_sequence_computation(const HiddenSemiMarkov &input,
      const MarkovianSequences &seq, bool characteristic_flag = true)
  {
    std:stringstream os;
    SIMPLE_METHOD_TEMPLATE_1(input, state_sequence_computation,
        SemiMarkovData, os, seq, characteristic_flag);
  }

  static HiddenSemiMarkov*
  thresholding(const HiddenSemiMarkov &input, double min_probability)
  {
    return input.thresholding(min_probability);
  }

  static SemiMarkovData*
  simulation_markovian_sequences(const HiddenSemiMarkov &input,
      int nb_sequence, const MarkovianSequences input_seq, bool counting_flag)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData, nb_sequence,
        input_seq, counting_flag);
  }

  static SemiMarkovData*
  simulation_histogram(const HiddenSemiMarkov &input,
      const FrequencyDistribution &hlength, bool counting_flag, bool divergence_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData, hlength,
        counting_flag, divergence_flag);

  }

  static SemiMarkovData*
  simulation_nb_sequences(const HiddenSemiMarkov &input, int nb_sequence,
      int length, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData, nb_sequence,
        length, counting_flag);
  }


  static DistanceMatrix*
  divergence_computation_histo(const HiddenSemiMarkov &input,
      boost::python::list input_markov,
      boost::python::list input_histogram_length, const char *filename)
  {
    HEADER_OS(DistanceMatrix);
    CREATE_ARRAY(input_markov, const HiddenSemiMarkov*, markov);
    CREATE_ARRAY(input_histogram_length, FrequencyDistribution*, hlength);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        hlength.get(), filename);
    FOOTER_OS;
  }

  static DistanceMatrix*
  divergence_computation_length(const HiddenSemiMarkov &input,
      boost::python::list input_markov, int nb_sequence, int length,
      const char *filename)
  {
    HEADER_OS(DistanceMatrix);
    CREATE_ARRAY(input_markov, const HiddenSemiMarkov*, markov);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_sequence, length, filename);
    FOOTER_OS;
  }


  static DistanceMatrix*
  divergence_computation_sequences(const HiddenSemiMarkov &input,
      boost::python::list &input_markov, boost::python::list &input_sequence,
      int nb_seq, char *filename)
  {
    // there is as much VariableOrderMarkov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    HEADER_OS(DistanceMatrix);

    CREATE_ARRAY(input_markov, const HiddenSemiMarkov*, markov);
    CREATE_ARRAY(input_sequence, const MarkovianSequences*, sequence);

    ret = input.divergence_computation(error, os, markov_size + 1, markov.get(),
        sequence_size, sequence.get(), filename);

    FOOTER_OS;
  }

/*  static MultiPlotSet*
  get_plotable(const HiddenSemiMarkov &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }
*/
};



void class_hidden_semi_markov() {



  class_<HiddenSemiMarkov, bases<SemiMarkov > >
    ("_HiddenSemiMarkov", "HiddenSemiMarkov", no_init)
    .def(init <HiddenSemiMarkov, bool, int>())
    .def("__init__", make_constructor(WRAP::hidden_semi_markov_from_file))

    .def(self_ns::str(self)) //__str__
    .def("file_ascii_write", WRAP::file_ascii_write, "Save vector summary into a file")

    DEF_RETURN_VALUE("state_sequence_computation", WRAP::state_sequence_computation, args(""),"")
    DEF_RETURN_VALUE("thresholding", WRAP::thresholding, args("index"), "todo")
    DEF_RETURN_VALUE("simulation_histogram", &WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", &WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_markovian_sequences", &WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")

    DEF_RETURN_VALUE("divergence_computation_histo", WRAP::divergence_computation_histo, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_length", WRAP::divergence_computation_length, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_sequences", WRAP::divergence_computation_sequences, args("input", "input_markov", "input_sequence", "filename"), "todo")
  //  DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")

;

/*
 *
    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,             int inb_output_process , NonparametricProcess **pobservation ,    int length , bool counting_flag)  :SemiMarkov(pchain , poccupancy , inb_output_process , pobservation , length ,  counting_flag) {}
    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,          int inb_output_process , NonparametricProcess **nonparametric_observation ,    DiscreteParametricProcess **parametric_observation , int length , bool counting_flag)    :SemiMarkov(pchain , poccupancy , inb_output_process , nonparametric_observation ,     parametric_observation , length , counting_flag) {}
    HiddenSemiMarkov(const HiddenSemiMarkov &hsmarkov , bool data_flag = true ,     int param = I_DEFAULT)

    bool spreadsheet_write(StatError &error , const char *path) const;

    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,  int inb_output_process , NonparametricProcess **pobservation ,  int length , bool counting_flag)    :SemiMarkov(pchain , poccupancy , inb_output_process , pobservation , length ,  counting_flag) {}
    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,    int inb_output_process , NonparametricProcess **nonparametric_observation ,int observationn , int length , bool counting_flag)    :SemiMarkov(pchain , poccupancy , inb_output_process , nonparametric_observation ,                 parametric_observation , length , counting_flag) {}


    double likelihood_computation(const MarkovianSequences &seq , double *posterior_probability = 0 ,  int index = I_DEFAULT) const;

    bool state_profile_write(StatError &error , std::ostream &os , const MarkovianSequences &iseq ,   int identifier = I_DEFAULT , int output = SSTATE ,  char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,  int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const char *path , const MarkovianSequences &iseq ,  int identifier = I_DEFAULT , int output = SSTATE ,  char format = 'a' , int state_sequence = GENERALIZED_VITERBI , int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_ascii_write(StatError &error , std::ostream &os , int identifier ,int output = SSTATE , int state_sequence = GENERALIZED_VITERBI ,int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const char *path ,  int identifier = I_DEFAULT , int output = SSTATE , char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,  int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_plot_write(StatError &error , const char *prefix , const MarkovianSequences &iseq , int identifier ,int output = SSTATE , const char *title = 0) const;
    bool state_profile_plot_write(StatError &error , const char *prefix , int identifier ,   int output = SSTATE , const char *title = 0) const;

    SemiMarkovData* state_sequence_computation(StatError &error , ostream &os , const MarkovianSequences &seq , bool characteristic_flag = true) const;
*/
}


#undef WRAP
