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
#include "sequence_analysis/semi_markov.h"
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



#define WRAP SemiMarkovWrap
class SemiMarkovWrap
{

public:

  static boost::shared_ptr<SemiMarkov>
  semi_markov_from_file(char* filename, int length, bool counting_flag, double cumul_threshold)
  {
    StatError error;

    SemiMarkov *semi_markov = NULL;

    semi_markov = semi_markov_ascii_read(error, filename, length,
        counting_flag, cumul_threshold);


    return boost::shared_ptr<SemiMarkov>(semi_markov);
  }

  static SemiMarkovData*
  extract_data(const SemiMarkov &input)
  {
    StatError err;
    SemiMarkovData *ret;
    ret = input.extract_data(err);
    return ret;
    //SIMPLE_METHOD_TEMPLATE_1(input, extract_data, SemiMarkovData);
  }

  static DiscreteParametricModel*
  extract_histogram(const SemiMarkov &input, int state, int histogram_type)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteParametricModel, state, histogram_type);
  }

  static DiscreteParametricModel*
  extract(const SemiMarkov &input, int type, int variable, int value)
  {

    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteParametricModel, type, variable,
        value);
  }



  static DiscreteParametricModel*
  get_forward(const SemiMarkov &input)
  {
    // todo: output must be list of Forward
    // todo: export Forward in stat_tool!
    DiscreteParametricModel *ret = NULL;
    return ret;
  }

/*  static Forward*
  get_forward(const SemiMarkov &input, int state)
  {
    Forward *ret;
    ret = input.get_forward(state);
    return ret;
  }
  BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(get_forward_overloads, SemiMarkovWrap::get_forward, 0, 1); */


  static void
  file_ascii_write(const SemiMarkov &d, const char *path, bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }

  static void
  write_hidden_semi_markov_init_file(const SemiMarkov &d, const char *path)
  {
    bool result = true;
    StatError error;
    result = d.write_hidden_semi_markov_init_file(error, path);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }

  static SemiMarkov*
  thresholding(const SemiMarkov &input, double min_probability)
  {
	return input.thresholding(min_probability);
  }

  static SemiMarkovData*
  simulation_histogram(const SemiMarkov &input,
	   const FrequencyDistribution &hlength , bool counting_flag, bool divergence_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData,
   	    hlength, counting_flag, divergence_flag);
  }

  static SemiMarkovData*
  simulation_nb_elements(const SemiMarkov &input,
		 int nb_sequence , int length , bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData,
   	    nb_sequence, length, counting_flag);

  }

  static SemiMarkovData*
  simulation_markovian_sequences(const SemiMarkov &input,
		  int nb_sequence , const MarkovianSequences &iseq , bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, SemiMarkovData,
   	    nb_sequence, iseq, counting_flag);
  }


  static DistanceMatrix*
  divergence_computation_histo(const SemiMarkov &input,
      boost::python::list input_markov,
      boost::python::list input_histogram_length,  const char *filename)
  {
    HEADER_OS(DistanceMatrix);
    CREATE_ARRAY(input_markov, const SemiMarkov*, markov);
    CREATE_ARRAY(input_histogram_length, FrequencyDistribution*, hlength);
    ret = input.divergence_computation(error, os, markov_size, markov.get(), hlength.get(), filename);
    FOOTER_OS;
  }

  static DistanceMatrix*
  divergence_computation_length(const SemiMarkov &input,
     boost::python::list input_markov , int nb_sequence ,
     int length , const char *filename)
  {
    HEADER_OS(DistanceMatrix);
    CREATE_ARRAY(input_markov, const SemiMarkov*, markov);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_sequence, length, filename);
    FOOTER_OS;
  }


  static DistanceMatrix*
  divergence_computation_sequences(const SemiMarkov &input,
      boost::python::list &input_markov, boost::python::list &input_sequences,
      int nb_seq, const char *filename)
  {
    // there is as much VariableOrderMarkov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    HEADER_OS(DistanceMatrix);

    CREATE_ARRAY(input_markov, const SemiMarkov*, markov);
    CREATE_ARRAY(input_sequences, const MarkovianSequences*, sequences);

    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_seq, sequences.get(), filename);
    FOOTER_OS;
  }

  static MultiPlotSet*
  get_plotable(const SemiMarkov &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }



};

// Boost declaration

void
class_semi_markov()
{

  class_<SemiMarkov, bases<StatInterface> >
  ("_SemiMarkov",  "SemiMarkov\n"
      "Constructors from a file required 3 arguments: length(int), counting_flag(boolean) and cumul_threshold (double)")
    .def("__init__", make_constructor(SemiMarkovWrap::semi_markov_from_file))

    .def(self_ns::str(self)) //__str__

    .add_property("nb_iterator", &SemiMarkov::get_nb_iterator, "returns nb iterator")
    .add_property("nb_output_process", &SemiMarkov::get_nb_output_process, "returns nb output process")

    //DEF_RETURN_VALUE("get_parametric_process", &SemiMarkov::get_parametric_process, args("variable_index"), "returns parametric process corresponding to the given variable")
    .def("extract_histogram", SemiMarkovWrap::extract_histogram, return_value_policy< manage_new_object >(), "todo")
    .def("extract", SemiMarkovWrap::extract, return_value_policy< manage_new_object >(), "todo")

    DEF_RETURN_VALUE_NO_ARGS("get_semi_markov_data", &SemiMarkov::get_semi_markov_data, "returns semi_markov_data")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", SemiMarkovWrap::extract_data, "returns semi_markov_data")

    .def("file_ascii_write", SemiMarkovWrap::file_ascii_write,"Save vector summary into a file")
    .def("write_hidden_semi_markov_init_file", SemiMarkovWrap::write_hidden_semi_markov_init_file, "save a SMM to a file in order to create a initialization file for HSMM")
    DEF_RETURN_VALUE("thresholding", SemiMarkovWrap::thresholding, args("index"), "todo")

    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_nb_elements",WRAP::simulation_nb_elements, args("todo"), "simulation")
    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("todo"), "simulation")

    DEF_RETURN_VALUE("divergence_computation_histo", WRAP::divergence_computation_histo, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_length", WRAP::divergence_computation_length, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_sequences", WRAP::divergence_computation_sequences, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")

    ;


  //todo file_ascii_write not accessible ?
  /*
   *
   SemiMarkov();
   SemiMarkov(char itype , int inb_state , int inb_output_process , int *nb_value);
   SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy , const NonparametricProcess *pobservation , int length ,  bool counting_flag);
   SemiMarkov(const SemiMarkov &smarkov , bool data_flag = true ,  int param = I_DEFAULT)  :Chain(smarkov) { copy(smarkov , data_flag , param); }

    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,const char *title = 0) const;

   void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
   void characteristic_computation(const SemiMarkovData &seq , bool counting_flag, int variable = I_DEFAULT , bool length_flag = true);

   double likelihood_computation(const MarkovianSequences &seq , int index) const;
   double likelihood_computation(const SemiMarkovData &seq) const;


   Forward** get_forward() const { return forward; }

   NonparametricSequenceProcess* get_nonparametric_process(int variable)   const { return nonparametric_process[variable]; }
   DiscreteParametricProcess** get_parametric_process() const { return parametric_process; }


   */

  ;
}
#undef WRAP


class SemiMarkovDataWrap
{

public:

   static SemiMarkovData*
   remove_index_parameter(const SemiMarkovData &input)
   {
     SIMPLE_METHOD_TEMPLATE_0(input, remove_index_parameter, SemiMarkovData);
   }

   static DiscreteDistributionData*
   extract(const SemiMarkovData &input, int type, int variable, int value)
   {
     SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteDistributionData, type, variable, value);
   }

 static MarkovianSequences*
 build_auxiliary_variable(const SemiMarkovData &input)
  {
    StatError error;
    MarkovianSequences* ret;
    ret = input.build_auxiliary_variable(error);
    if (!ret) 
      sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }


};

void
class_semi_markov_data()
{

  class_<SemiMarkovData, bases<MarkovianSequences> > ("_SemiMarkovData", "SemiMarkovData")
    .def(init<MarkovianSequences> ())
    .def(init< MarkovianSequences, char, bool> ())
    .def(init<SemiMarkovData, bool, char> ())

    .add_property("likelihood", &SemiMarkovData::get_likelihood, "returns likelihood")
    .add_property("restoration_likelihood", &SemiMarkovData::get_restoration_likelihood, "returns restoration likelihood")

    .def("get_posterior_probability", &SemiMarkovData::get_posterior_probability, args("index"))

    DEF_RETURN_VALUE_NO_ARGS("get_semi_markov", &SemiMarkovData::get_semi_markov, "returns semi_markov")
    DEF_RETURN_VALUE_NO_ARGS("get_chain_data", &SemiMarkovData::get_chain_data, "returns chain data")
    DEF_RETURN_VALUE_NO_ARGS("remove_index_parameter", &SemiMarkovDataWrap::remove_index_parameter, "remove index parameter")
    DEF_RETURN_VALUE_NO_ARGS("build_auxiliary_variable", &SemiMarkovDataWrap::build_auxiliary_variable, "calls build_auxialiary_varibles methods and returns markovian sequences object")
    

 ;

  // DONE
  // todo do we need the constructor and build_transition?
  /*
   SemiMarkovData();
   SemiMarkovData(const FrequencyDistribution &ihlength , int inb_variable , bool init_flag = false);
   void build_transition_count(const SemiMarkov *smarkov = 0);
   */

}

class SemiMarkovIteratorWrap
{

public:



  static boost::python::list
  simulation(SemiMarkovIterator &input, int nb_sequence = 1,
      bool initialisation=false)
  {
    StatError error;
    int **sequence;

    SemiMarkov *sm;
    sm = input.get_semi_markov();

    sequence = input.simulation(nb_sequence, initialisation);
    boost::python::list output_sequence;

    int i, j;
    for (j=0; j < sm->get_nb_output_process()+1; j++)
    {
        boost::python::list line;
        for (i=0; i < nb_sequence; i++)
        {
            line.append(sequence[j][i]);
        }
        output_sequence.append(line);
    }

    // delete the sequence allocated in C++
    for (j=0; j < sm->get_nb_output_process()+1; j++)
        {
        delete [] sequence[j];
        sequence[j] = 0;
        }
    delete [] sequence;

    return output_sequence;
  }
};


void
class_semi_markov_iterator()
{

  class_<SemiMarkovIterator, boost::shared_ptr<SemiMarkovIterator> >
  ("_SemiMarkovIterator", "SemiMarkovIterator", no_init)

  .def(init<SemiMarkov*>() [with_custodian_and_ward_postcall<1, 2>()])


  .add_property("state", &SemiMarkovIterator::get_state)
  .add_property("occupancy", &SemiMarkovIterator::get_occupancy)
  .add_property("counter", &SemiMarkovIterator::get_counter)
  .add_property("nb_variable", &SemiMarkovIterator::get_nb_variable)


  .def("simulation", SemiMarkovIteratorWrap::simulation,
        args("nb_sequence","initialisation"), "simulation")


  ;
}
