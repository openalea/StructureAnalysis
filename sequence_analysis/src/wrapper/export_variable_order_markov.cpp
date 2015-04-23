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



#define WRAP VariableOrderMarkovWrap
class VariableOrderMarkovWrap {

public:

  static boost::shared_ptr<VariableOrderMarkov>
  variable_order_markov_from_file(char *filename, int length)
  {
    StatError error;
    VariableOrderMarkov *vom = NULL;
    vom = variable_order_markov_ascii_read(error, filename, length);
    if (!vom)
      {
        sequence_analysis::wrap_util::throw_error(error);
      }
    return boost::shared_ptr<VariableOrderMarkov>(vom);
  }


  static void
  file_ascii_write(const VariableOrderMarkov &d, const char *path,
      bool exhaustive)
  {
    bool result = true;
    StatError error;

    result = d.ascii_write(error, path, exhaustive);
    if (!result)
      sequence_analysis::wrap_util::throw_error(error);
  }

  static DiscreteParametricModel*
  extract(const VariableOrderMarkov &input, int type, int variable, int value)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteParametricModel, type, variable,
        value);
  }

  static VariableOrderMarkovData*
  extract_data(const VariableOrderMarkov &input)
  {
    SIMPLE_METHOD_TEMPLATE_0(input, extract_data, VariableOrderMarkovData);
  }

  static Correlation*
  output_autocorrelation_computation(const VariableOrderMarkov &input,
      int variable, int output, int max_lag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, output_autocorrelation_computation,
        Correlation, variable, output, max_lag);
  }

  static Correlation*
  state_autocorrelation_computation(const VariableOrderMarkov &input,
      int state, int max_lag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, state_autocorrelation_computation,
        Correlation, state, max_lag);
  }

  static double
  likelihood_computation(const VariableOrderMarkov &input,
      const MarkovianSequences &ms, int index)
  {
    return input.likelihood_computation(ms, index);
  }

  static VariableOrderMarkov*
  thresholding(const VariableOrderMarkov &input, double min_probability)
  {
    return input.thresholding(min_probability);
  }

  static DistanceMatrix*
  divergence_computation_histo(const VariableOrderMarkov &input,
      boost::python::list input_markov,
      boost::python::list input_histogram_length, const char *filename)
  {
    HEADER_OS(DistanceMatrix);CREATE_ARRAY(input_markov, const VariableOrderMarkov*, markov);
    CREATE_ARRAY(input_histogram_length, FrequencyDistribution*, hlength);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        hlength.get(), filename);
    FOOTER_OS;
  }

  static DistanceMatrix*
  divergence_computation_length(const VariableOrderMarkov &input,
      boost::python::list input_markov, int nb_sequence, int length,
      const char *filename)
  {
    HEADER_OS(DistanceMatrix);CREATE_ARRAY(input_markov, const VariableOrderMarkov*, markov);
    ret = input.divergence_computation(error, os, markov_size, markov.get(),
        nb_sequence, length, filename);
    FOOTER_OS;
  }

  static DistanceMatrix*
  divergence_computation_sequences(const VariableOrderMarkov &input,
      boost::python::list &input_markov, boost::python::list &input_sequence,
      int nb_seq, char *filename)
  {
    // there is as much VariableOrderMarkov elmts as Markovian elts
    // 1 for input and N-1 for input_markov and N input_sequence

    HEADER_OS(DistanceMatrix);

    CREATE_ARRAY(input_markov, const VariableOrderMarkov*, markov);
    CREATE_ARRAY(input_sequence, const MarkovianSequences*, sequence);

    ret = input.divergence_computation(error, os, markov_size + 1, markov.get(),
        sequence_size, sequence.get(), filename);

    FOOTER_OS;
  }

  static VariableOrderMarkovData*
  simulation_markovian_sequences(const VariableOrderMarkov &input,
      int nb_sequence, const MarkovianSequences input_seq, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, VariableOrderMarkovData,
        nb_sequence, input_seq, counting_flag);
  }

  static VariableOrderMarkovData*
  simulation_histogram(const VariableOrderMarkov &input,
      const FrequencyDistribution &hlength, bool counting_flag, bool divergence_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, VariableOrderMarkovData,
        hlength, counting_flag, divergence_flag);
  }

  static VariableOrderMarkovData*
  simulation_nb_sequences(const VariableOrderMarkov &input, int nb_sequence,
      int length, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, VariableOrderMarkovData,
        nb_sequence, length, counting_flag);
  }

  static MultiPlotSet*
  get_plotable(const VariableOrderMarkov &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
    }


};

// Boost declaration

void class_variable_order_markov() {

  class_<VariableOrderMarkov, bases<StatInterface > >
  ("_VariableOrderMarkov", "The only possible constructor is VariableOrderMarkov(filename) \n"
      "The following constructors will be available once the Chain class is made public (protected right now):\n"
      " * constructor(type(char in ['o','e']), nb_state(int), nb_row(int)\n"
      " * constructor(type(char in ['o','e']), nb_state(int), nb_row(int), max_order(int)\n"
      " * constructor(type(char in ['o','e']), nb_state(int), order(int), init_flag(bool), output_process=0, nb_value=0\n", no_init
  )
    .def("__init__", make_constructor(VariableOrderMarkovWrap::variable_order_markov_from_file))

      // type = 'o' : ordinaire, or 'e' : en equilibre des probabilites de transition
    //.def(init <char, int, int>())
    //.def(init <char, int, int,int>())
    //.def(init <char, int, int, bool, optional<int, int> >())

    .def(self_ns::str(self)) //__str__

    .def("file_ascii_write", WRAP::file_ascii_write,"Save vector summary into a file")
//    .def("ascii_data_write", WRAP::ascii_data_write,"Return a string with the object representation")
//    .def("file_ascii_data_write", WRAP::file_ascii_data_write,"Save vector data into a file")

/*    .add_property("nb_iterator", &VariableOrderMarkov::get_nb_iterator, "todo")
    .add_property("max_order", &VariableOrderMarkov::get_max_order, "todo")
    .add_property("nb_output_process", &VariableOrderMarkov::get_nb_output_process, "todo")

    .def("get_memory_type", &VariableOrderMarkov::get_memory_type, args("memory"),"todo")
    .def("get_order", &VariableOrderMarkov::get_order, args("memory"),"todo")
    .def("get_state", &VariableOrderMarkov::get_state, args("memory", "lag"),"todo")
    .def("get_parent", &VariableOrderMarkov::get_parent, args("memory"),"todo")
    .def("get_child", &VariableOrderMarkov::get_child, args("memory", "state"),"todo")
    .def("get_next", &VariableOrderMarkov::get_next, args("memory", "state"),"todo")
    .def("get_nb_memory", &VariableOrderMarkov::get_nb_memory, args("memory"),"todo")
    .def("get_previous", &VariableOrderMarkov::get_previous,args("memory", "state"), "todo") */

    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data, "returns variable_order_markov_data")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")

    DEF_RETURN_VALUE("extract", WRAP::extract, args("type","variable","value"), "returns parametric model")
    DEF_RETURN_VALUE_NO_ARGS("extract_data", WRAP::extract_data, "returns variable_order_markov_data")
    DEF_RETURN_VALUE("output_autocorrelation_computation", WRAP::output_autocorrelation_computation, args("variable", "output", "max_lag"), "todo")
    DEF_RETURN_VALUE("state_autocorrelation_computation", WRAP::state_autocorrelation_computation, args("state", "max_lag"), "todo")

    .def("likelihood_computation", WRAP::likelihood_computation, args("seq", "index"),"todo")

    DEF_RETURN_VALUE("thresholding", WRAP::thresholding, args("index"), "todo")

    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")


    DEF_RETURN_VALUE("divergence_computation_histo", WRAP::divergence_computation_histo, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_length", WRAP::divergence_computation_length, args("input", "input_markov", "input_sequence", "filename"), "todo")
    DEF_RETURN_VALUE("divergence_computation_sequences", WRAP::divergence_computation_sequences, args("input", "input_markov", "input_sequence", "filename"), "todo")

    ;


/*
 *
 *
  VariableOrderMarkov(const VariableOrderMarkov &markov , int inb_output_process , int nb_value);
  VariableOrderMarkov(const VariableOrderMarkov *pmarkov ,  const NonparametricProcess *pobservation , int length);
  VariableOrderMarkov(const VariableOrderMarkov &markov , bool data_flag = true)  :Chain(markov) { copy(markov , data_flag); }

  void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
  void characteristic_computation(const VariableOrderMarkovData &seq , bool counting_flag ,
  int variable = I_DEFAULT , bool length_flag = true);

  double likelihood_computation(const VariableOrderMarkovData &seq) const;


  VariableOrderMarkovData* get_markov_data() const { return markov_data; }
  NonparametricSequenceProcess* get_nonparametric_process(int variable) const{ return nonparametric_process[variable]; }
  DiscreteParametricProcess** get_parametric_process() const { return parametric_process; }
  DiscreteParametricProcess* get_parametric_process(int variable)const { return parametric_process[variable]; }
*/




}
#undef WRAP


#define WRAP VariableOrderMarkovDataWrap
class VariableOrderMarkovDataWrap {

public:

 static MarkovianSequences*
 build_auxiliary_variable(const VariableOrderMarkovData &input)
  {
    StatError error;
    MarkovianSequences* ret;
    ret = input.build_auxiliary_variable(error);
    if (!ret) 
      sequence_analysis::wrap_util::throw_error(error);
    return ret;
  }


};

void class_variable_order_markov_data() {

    class_<VariableOrderMarkovData, bases<MarkovianSequences > >
    ("_VariableOrderMarkovData", "VariableOrderMarkovData")

    DEF_RETURN_VALUE_NO_ARGS("build_auxiliary_variable", WRAP::build_auxiliary_variable, 
       "calls build_auxialiary_varibles methods and returns markovian sequences object")
    ;


        /*
   VariableOrderMarkovData();
   VariableOrderMarkovData(const FrequencyDistribution &ihlength , int inb_variable , bool init_flag = false);
   VariableOrderMarkovData(const MarkovianSequences &seq);
   VariableOrderMarkovData(const MarkovianSequences &seq , char transform ,
                           bool initial_run_flag);
   VariableOrderMarkovData(const VariableOrderMarkovData &seq , bool model_flag = true ,
                           char transform = 'c')
   :MarkovianSequences(seq , transform) { copy(seq , model_flag); }
   ~VariableOrderMarkovData();
   VariableOrderMarkovData& operator=(const VariableOrderMarkovData &seq);

   DiscreteDistributionData* extract(StatError &error , int type ,
                                     int variable , int value) const;
   VariableOrderMarkovData* remove_index_parameter(StatError &error) const;

   Correlation* state_autocorrelation_computation(StatError &error , int istate ,
                                                  int max_lag = MAX_LAG) const;
   Correlation* output_autocorrelation_computation(StatError &error , int variable ,
                                                   int output , int max_lag = MAX_LAG) const;

   std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                  bool exhaustive = false) const;
   bool ascii_data_write(StatError &error , const char *path ,
                         char format = 'c' , bool exhaustive = false) const;

   std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
   bool ascii_write(StatError &error , const char *path ,
                    bool exhaustive = false) const;
   bool spreadsheet_write(StatError &error , const char *path) const;
   bool plot_write(StatError &error , const char *prefix ,
                   const char *title = 0) const;

   void build_transition_count(const VariableOrderMarkov &markov ,
                               bool begin = true , bool non_terminal = false);
   void order0_estimation(VariableOrderMarkov &markov) const;

   // acces membres de la classe

  VariableOrderMarkov* get_markov() const { return markov; }
   VariableOrderChainData* get_chain_data() const { return chain_data; }
   double get_likelihood() const { return likelihood; }
   double get_restoration_likelihood() const { return restoration_likelihood; }
   double get_posterior_probability(int index) const { return posterior_probability[index]; }
   */
}

#undef WRAP


class VariableOrderMarkovIteratorWrap
{

public:

  static boost::python::list
  simulation(VariableOrderMarkovIterator& input, int nb_sequence=1,
      bool initialisation=false)
  {
    StatError error;
    int **sequence;

    VariableOrderMarkov *sm;
    sm = input.get_markov();

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
class_variable_order_markov_iterator()
{

  class_<VariableOrderMarkovIterator >
  ("_VariableOrderMarkovIterator", "VariableOrderMarkovIterator", no_init)

  .def(init<VariableOrderMarkov*>()[with_custodian_and_ward_postcall<1, 2>()])

  .add_property("nb_variable", &VariableOrderMarkovIterator::get_nb_variable)
  .add_property("memory", &VariableOrderMarkovIterator::get_memory)

  .def("simulation", VariableOrderMarkovIteratorWrap::simulation,
      args("nb_sequence","initialisation"), "simulation")
    ;
}
