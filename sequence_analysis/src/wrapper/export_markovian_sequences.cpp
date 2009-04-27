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
#include "stat_tool/vectors.h"
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

using namespace boost::python;
using namespace boost;
using namespace stat_tool;




class MarkovianSequencesWrap {

public:

   static Markovian_sequences* markovian_sequences1(
            int nb_sequence,
            boost::python::list& input_identifiers,
            boost::python::list& input_ilength,
            int input_index_parameter_type,
            int input_nb_variable,
            boost::python::list& input_itype,
            bool init_flag)
   {
    Markovian_sequences *seq = NULL;
    // todo
    return seq;
    }

    static Markovian_sequences* markovian_sequences2(
            int nb_sequence,
            boost::python::list& input_identifiers,
            boost::python::list& input_ilength,
            int input_index_parameter_type,
            int input_nb_variable,
            boost::python::list& input_itype,
            bool init_flag)
   {
    Markovian_sequences *seq = NULL;
    // todo
    return seq;
    }

    static Markovian_sequences* markovian_sequences3(
            int nb_sequence,
            boost::python::list& input_identifiers,
            boost::python::list& input_ilength,
            int input_index_parameter_type,
            int input_nb_variable,
            boost::python::list& input_itype,
            bool init_flag)
   {
    Markovian_sequences *seq = NULL;
    // todo
    return seq;
    }


    static Distribution_data* extract(const Markovian_sequences& seq, int type, int variable, int value)
    {
        Format_error error;
        Distribution_data* ret;
        ret = seq.extract(error, type, variable, value);
        if (!ret)
            stat_tool::wrap_util::throw_error(error);
        return ret;
    }

    // Merge
  	static Markovian_sequences* merge(const Markovian_sequences& input_seq,
   			const boost::python::list& seqs) {
  		Format_error error;
   		Markovian_sequences * ret = NULL;

   		int nb_seq = len(seqs);
   		stat_tool::wrap_util::auto_ptr_array<const Markovian_sequences *> sequens(
   			new const Markovian_sequences*[nb_seq]);

   		for (int i = 0; i < nb_seq; i++)
   			sequens[i] = boost::python::extract< Markovian_sequences *> (seqs[i]);

   		ret = input_seq.merge(error, nb_seq, sequens.get());
   		if (!ret)
   			stat_tool::wrap_util::throw_error(error);

   		return ret;
   	}

  	static Markovian_sequences* cluster(const Markovian_sequences& seq,
  			int variable, int step, int mode = FLOOR)
  	{
  	    Format_error error;
  	    Markovian_sequences* ret;
  	    ret = seq.cluster(error, variable, step, mode);
  	    if (!ret)
  	        stat_tool::wrap_util::throw_error(error);
  	    return ret;
  	}


  	static Markovian_sequences* transcode(const Markovian_sequences& seq,
		int variable, boost::python::list& symbol, bool add_flag = false) {

  		Format_error error;

  		int nb_symbol = len(symbol);
  		stat_tool::wrap_util::auto_ptr_array<int> l(new int[nb_symbol]);

  		//? needed ?
  		int expected_nb_symbol = (int) (seq.get_max_value(variable - 1)
  				- seq.get_min_value(variable - 1)) + 1;
  		//? needed ?
  		if (nb_symbol != expected_nb_symbol)
  			stat_tool::wrap_util::throw_error("Bad number of Symbol");

  		for (int i = 0; i < nb_symbol; i++)
  			l[i] = boost::python::	extract<int> (symbol[i]);

  		Markovian_sequences* ret = seq.transcode(error, variable, l.get(), add_flag);

  		if (!ret)
  			stat_tool::wrap_util::throw_error(error);

  		return ret;
	}

	static Markovian_sequences* remove_index_parameter(const Markovian_sequences& seq)
  	{
  	    Format_error error;
  	    Markovian_sequences* ret;
  	    ret = seq.remove_index_parameter(error);
  	    if (!ret)
  	        stat_tool::wrap_util::throw_error(error);
  	    return ret;
  	}




};

// Boost declaration

void class_markovian_sequences() {

	class_<Markovian_sequences, bases<Sequences> > ("_Markovian_sequences", "Markovian_sequences")
    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences1))
    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences2))
    .def("__init__", make_constructor(MarkovianSequencesWrap::markovian_sequences3))
    .def(init <Sequences>())
    .def(init <Markovian_sequences,char, int>())
	.def("extract", MarkovianSequencesWrap::extract,
		return_value_policy<manage_new_object> (),
        python::args("type", "variable","value"),
	   "Extract distribution data")
	.def("merge", &MarkovianSequencesWrap::merge,
		return_value_policy<manage_new_object> (),
		python::args("sequences"),
	    "Merge")
	.def("cluster", &MarkovianSequencesWrap::cluster,
		return_value_policy<manage_new_object> (),
		python::args("variable", "step", "mode"),
	    "Cluster")
	.def("transcode", &MarkovianSequencesWrap::transcode,
		return_value_policy<manage_new_object> (),
		python::args("variable", "symbol", "add_flag"),
	    "Transcode")
	 .def("remove_index_parameter", &MarkovianSequencesWrap::remove_index_parameter,
	    return_value_policy<manage_new_object> (),
	    "Remove index parameter")

/*
    Markovian_sequences* consecutive_values(Format_error &error , std::ostream &os ,
                                            int ivariable , bool add_flag = false) const;
    Markovian_sequences* cluster(Format_error &error , int ivariable , int nb_class ,
                                 int *ilimit , bool add_flag = false) const;
    Markovian_sequences* cluster(Format_error &error , int variable , int nb_class ,
                                 double *ilimit) const;

    Markovian_sequences* select_variable(Format_error &error , int inb_variable ,
                                         int *ivariable , bool keep = true) const;
    Markovian_sequences* merge_variable(Format_error &error , int nb_sample ,
                                        const Markovian_sequences **iseq , int ref_sample = I_DEFAULT) const;
    Markovian_sequences* remove_variable_1() const;

    Markovian_sequences* initial_run_computation(Format_error &error) const;
    Markovian_sequences* add_absorbing_run(Format_error &error ,
                                           int sequence_length = I_DEFAULT ,
                                           int run_length = I_DEFAULT) const;

    Markovian_sequences* split(Format_error &error , int step) const;
i Markovian_sequences* split(Format_error &error , int step) const;

    bool ascii_data_write(Format_error &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

    bool transition_count(Format_error &error , std::ostream &os , int max_order ,
                          bool begin = false , int estimator = MAXIMUM_LIKELIHOOD ,
                          const char *path = 0) const;
    bool word_count(Format_error &error , std::ostream &os , int variable , int word_length ,
                    int begin_state = I_DEFAULT , int end_state = I_DEFAULT ,
                    int min_frequency = 1) const;
    bool mtg_write(Format_error &error , const char *path , int *itype) const;


*/

;
}




void class_self_transition() {

	class_<Self_transition, bases<Curves> > ("_Self_transition", "Self_transition", no_init)
	.def(init <int>())
	.def("monomolecular_regression", &Self_transition::monomolecular_regression,
	    return_value_policy<manage_new_object> (),
	   "return function monomolecular regression")
	.def("logistic regression", &Self_transition::logistic_regression,
		return_value_policy<manage_new_object> (),
		"return function logistic regression")
    ;

}
