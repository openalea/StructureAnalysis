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
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/distribution.h"
#include "stat_tool/markovian.h"
#include "stat_tool/vectors.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/stat_label.h"

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
using namespace stat_tool;
using namespace sequence_analysis;



#define WRAP NonHomogeneousMarkovWrap
class NonHomogeneousMarkovWrap
{

public:

  static boost::shared_ptr<NonhomogeneousMarkov>
  nonhomogeneous_markov_from_file(char *filename, int length)
  {
    //olf_format should be true
    StatError error;
    NonhomogeneousMarkov *nonhomo = NULL;
    nonhomo = nonhomogeneous_markov_ascii_read(error, filename, length);
    /*if (!nonhomo)
      {
        sequence_analysis::wrap_util::throw_error(error);
      }
    */
    return boost::shared_ptr<NonhomogeneousMarkov>(nonhomo);
  }

  static NonhomogeneousMarkovData*
  simulation_markovian_sequences(const NonhomogeneousMarkov &input,
      int nb_sequence, const MarkovianSequences input_seq, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, NonhomogeneousMarkovData,
        nb_sequence, input_seq, counting_flag);
  }

  static NonhomogeneousMarkovData*
  simulation_histogram(const NonhomogeneousMarkov &input,
      const FrequencyDistribution &hlength, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, NonhomogeneousMarkovData,
        hlength, counting_flag);
  }

  static NonhomogeneousMarkovData*
  simulation_nb_sequences(const NonhomogeneousMarkov &input, int nb_sequence,
      int length, bool counting_flag)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, simulation, NonhomogeneousMarkovData,
        nb_sequence, length, counting_flag);
  }

  static DiscreteParametricModel*
  extract(const NonhomogeneousMarkov &input, int type, int state)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract, DiscreteParametricModel,
        type, state);
  }

 static MultiPlotSet*
  get_plotable(const NonhomogeneousMarkov &p)
  {
    StatError error;
    MultiPlotSet* ret = p.get_plotable();
    if (!ret)
      ERROR;
    return ret;
  }

};

// Boost declaration

void class_nonhomogeneous_markov() {

  class_<NonhomogeneousMarkov, bases<StatInterface> > ("_NonHomogeneousMarkov", "NonHomogeneousMarkov")
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_nb_state_and_ident_list))
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_chain_and_self_transition))
    .def("__init__", make_constructor(NonHomogeneousMarkovWrap::nonhomogeneous_markov_from_file))

    .def(self_ns::str(self)) //__str__

    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")
    .def("get_homogeneity", &NonhomogeneousMarkov::get_homogeneity, python::args("index"),"return homogeneity")
    .def("get_self_transition", &NonhomogeneousMarkov::get_self_transition, return_value_policy<manage_new_object> (),	python::args("index"),"return Function* of self transition")

    DEF_RETURN_VALUE("simulation_histogram", WRAP::simulation_histogram, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_nb_sequences", WRAP::simulation_nb_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE("simulation_markovian_sequences", WRAP::simulation_markovian_sequences, args("nb_sequence", "input_seq", "counting_flag"), "todo")
    DEF_RETURN_VALUE_NO_ARGS("get_plotable", WRAP::get_plotable, "Return a plotable")


    ;

	//.def("get_process", &NonhomogeneousMarkov::get_process,
	//	"return non parametric sequence process")

/*
    NonhomogeneousMarkov();
    NonhomogeneousMarkov(int inb_state , int *ident);
    NonhomogeneousMarkov(const Chain *pchain , const Function **pself_transition , int length);
    NonhomogeneousMarkov(const NonhomogeneousMarkov &markov , bool data_flag = true , bool haracteristic_flag = true)    :Chain(markov) { copy(markov , data_flag , characteristic_flag); }


    void characteristic_computation(int length , bool counting_flag);
    void characteristic_computation(const NonhomogeneousMarkovData &seq , bool counting_flag ,            bool length_flag = true);

    double likelihood_computation(const MarkovianSequences &seq , int index = I_DEFAULT) const;

    NonhomogeneousMarkovData* get_markov_data() const { return markov_data; }
    NonparametricSequenceProcess* get_process() const { return process; }
   */
}
#undef WRAP



#define WRAP NonHomogeneousMarkovDataWrap
class NonHomogeneousMarkovDataWrap
{

public:
  static DiscreteDistributionData*
  extract(const NonhomogeneousMarkovData &input, int type, int state)
  {
    SIMPLE_METHOD_TEMPLATE_1(input, extract,
        DiscreteDistributionData, type, state);
  }

  static NonhomogeneousMarkovData*
  remove_index_parameter(const NonhomogeneousMarkovData& input)
  {
    SIMPLE_METHOD_TEMPLATE_0(input, remove_index_parameter,
        NonhomogeneousMarkovData);
  }

};

void class_nonhomogeneous_markov_data()
{
  class_<NonhomogeneousMarkovData, bases<MarkovianSequences > >
  ("_NonHomogeneousMarkovData", "NonHomogeneousMarkovData")

    .def(init<MarkovianSequences& >())

    .add_property("likelihood", &NonhomogeneousMarkovData::get_likelihood, "returns likelihood")
    .def("extract", WRAP::extract, return_value_policy<manage_new_object> (),  python::args("type", "state"), "Extract distribution data")

  ;
  /*
    NonhomogeneousMarkovData();
    NonhomogeneousMarkovData(const FrequencyDistribution &ihlength);
    NonhomogeneousMarkovData(const NonhomogeneousMarkovData &seq , bool model_flag = true , char transform = 'c') :MarkovianSequences(seq , transform) { copy(seq , model_flag); }

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,const char *title = 0) const;

    void build_transition_count();

    NonhomogeneousMarkov* get_markov() const { return markov; }
    ChainData* get_chain_data() const { return chain_data; }

    */

}

#undef WRAP
