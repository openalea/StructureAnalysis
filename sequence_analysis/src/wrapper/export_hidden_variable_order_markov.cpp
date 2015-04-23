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

#include "tool/config.h"

#include "stat_tool/stat_tools.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
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



#define WRAP HiddenVariableOrderMarkovWrap
class HiddenVariableOrderMarkovWrap {

public:

  static boost::shared_ptr<HiddenVariableOrderMarkov>
  hidden_variable_order_markov_from_file(char *filename, int length,
      double cumul_threshold)
  {
    StatError error;

    HiddenVariableOrderMarkov *hvom = NULL;

    hvom = hidden_variable_order_markov_ascii_read(error, filename, length,
        cumul_threshold);

    return boost::shared_ptr<HiddenVariableOrderMarkov>(hvom);
  }

  static HiddenVariableOrderMarkov*
  thresholding(const HiddenVariableOrderMarkov &input,
      double min_probability = MIN_PROBABILITY)
  {
    //todo check this function
    HiddenVariableOrderMarkov *ret = NULL;
    ret = input.thresholding(min_probability);
    return ret;
  }

  static VariableOrderMarkovData*
  simulation_markovian_sequences(const HiddenVariableOrderMarkov &input,
      int nb_sequence, const MarkovianSequences input_seq, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation,
           VariableOrderMarkovData, nb_sequence, input_seq, counting_flag);
  }

  static VariableOrderMarkovData*
  state_sequence_computation(const HiddenVariableOrderMarkov &input,
      const MarkovianSequences &seq, bool characteristic_flag = true)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, state_sequence_computation,
        VariableOrderMarkovData, seq, characteristic_flag);
  }

  static VariableOrderMarkovData*
  simulation_histogram(const HiddenVariableOrderMarkov &input,
      const FrequencyDistribution &hlength, bool counting_flag, bool divergence_flag)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, simulation, VariableOrderMarkovData,
        hlength, counting_flag, divergence_flag);
  }

  static VariableOrderMarkovData*
  simulation_nb_sequences(const HiddenVariableOrderMarkov &input,
      int nb_sequence, int length, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, VariableOrderMarkovData,
      nb_sequence, length, counting_flag);

    return ret;
  }



  static DistanceMatrix*
  divergence_computation_histo(const HiddenVariableOrderMarkov &input,
      boost::python::list input_markov,
      boost::python::list input_histogram_length, const char *filename)
  {
    HEADER_OS(DistanceMatrix);CREATE_ARRAY(input_markov, const HiddenVariableOrderMarkov*, markov);
    CREATE_ARRAY(input_histogram_length, FrequencyDistribution*, hlength);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        hlength.get(), filename);
    FOOTER_OS;
  }

  static DistanceMatrix*
  divergence_computation_length(const HiddenVariableOrderMarkov &input,
      boost::python::list input_markov, int nb_sequence, int length,
      const char *filename)
  {
    HEADER_OS(DistanceMatrix);CREATE_ARRAY(input_markov, const HiddenVariableOrderMarkov*, markov);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_sequence, length, filename);
    FOOTER_OS;
  }


  static DistanceMatrix*
  divergence_computation_sequences(const HiddenVariableOrderMarkov &input,
      boost::python::list &input_markov, boost::python::list &input_sequence,
      int nb_seq, char *filename)
  {
    // there is as much VariableOrderMarkov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    HEADER_OS(DistanceMatrix);

    CREATE_ARRAY(input_markov, const HiddenVariableOrderMarkov*, markov);
    CREATE_ARRAY(input_sequence, const MarkovianSequences*, sequence);

    ret = input.divergence_computation(error, os, markov_size + 1, markov.get(),
        sequence_size, sequence.get(), filename);

    FOOTER_OS;
  }

  static MultiPlotSet*
  state_profile_plotable_write(const HiddenVariableOrderMarkov &p, int identifier)
  {
    StatError error;
    MultiPlotSet* ret = p.state_profile_plotable_write(error, identifier);
    if (!ret)
      ERROR;
    return ret;
  }
  static MultiPlotSet*
   state_profile_plotable_write2(const HiddenVariableOrderMarkov &input,
         const MarkovianSequences &iseq, int identifier)
   {
     StatError error;
     MultiPlotSet* ret = input.state_profile_plotable_write(error, iseq, identifier);
     if (!ret)
       ERROR;
     return ret;
   }



};



void class_hidden_variable_order_markov() {

  class_<HiddenVariableOrderMarkov, bases<VariableOrderMarkov > >
    ("_HiddenVariableOrderMarkov", "HiddenVariableOrderMarkov")
    .def("__init__", make_constructor(WRAP::hidden_variable_order_markov_from_file))

    .def(self_ns::str(self)) //__str__

    DEF_RETURN_VALUE("thresholding", WRAP::thresholding, args("probability"), "todo")
    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("state_sequence_computation", WRAP::state_sequence_computation, args(""),"")


    DEF_RETURN_VALUE("divergence_computation_histo", WRAP::divergence_computation_histo, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_length", WRAP::divergence_computation_length, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_sequences", WRAP::divergence_computation_sequences, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("state_profile_plotable_write", WRAP::state_profile_plotable_write, args("identifier", "output"), "Return a plotable")
    DEF_RETURN_VALUE("state_profile_plotable_write2", WRAP::state_profile_plotable_write2, args("sequences", "identifier", "output"), "Return a plotable")
;
}
      /*
        HiddenVariableOrderMarkov() {}
          HiddenVariableOrderMarkov(const VariableOrderMarkov *pmarkov , int inb_output_process ,
                                       NonparametricProcess **nonparametric_observation ,
                                       DiscreteParametricProcess **parametric_observation , int length)

          HiddenVariableOrderMarkov(const HiddenVariableOrderMarkov &hmarkov , bool data_flag = true)

          std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
          bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
          bool spreadsheet_write(StatError &error , const char *path) const;

		 double likelihood_computation(const MarkovianSequences &seq ,double *posterior_probability = 0 ,  int index = I_DEFAULT) const;

          bool state_profile_write(StatError &error , std::ostream &os , const MarkovianSequences &iseq ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,            int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(StatError &error , const char *path , const MarkovianSequences &iseq ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_ascii_write(StatError &error , std::ostream &os , int identifier ,int state_sequence = GENERALIZED_VITERBI ,int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(StatError &error , const char *path ,int identifier = I_DEFAULT , char format = 'a' ,int state_sequence = GENERALIZED_VITERBI ,    int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_plot_write(StatError &error , const char *prefix ,  const MarkovianSequences &iseq ,    int identifier , const char *title = 0) const;
          bool state_profile_plot_write(StatError &error , const char *prefix ,   int identifier , const char *title = 0) const;




*/



#undef WRAP
