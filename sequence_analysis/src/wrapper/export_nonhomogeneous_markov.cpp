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

using namespace boost::python;
using namespace boost;
using namespace stat_tool;




class NonHomogeneousMarkovWrap {

public:


	  static Parametric_model* extract(const Nonhomogeneous_markov& seq,
			  int type, int state)
	    {
	        Format_error error;
	        Parametric_model* ret;
	        ret = seq.extract(error, type, state);
	        if (!ret)
	            stat_tool::wrap_util::throw_error(error);
	        return ret;
	    }


};

// Boost declaration

void class_nonhomogeneous_markov() {

	class_<Nonhomogeneous_markov, bases<STAT_interface , Chain > > ("_Nonhomogeneous_markov", "Nonhomogeneous_markov")
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_nb_state_and_ident_list))
    //.def("__init__", make_constructor(NonHomogeneousMarkovWrap::constructor_from_chain_and_self_transition))


	.def("extract", &NonHomogeneousMarkovWrap::extract,
		return_value_policy<manage_new_object> (),
        python::args("type", "state"),
	   "Extract distribution data")
	.def("get_homogeneity", &Nonhomogeneous_markov::get_homogeneity,
		python::args("index"),"return homogeneity")
	.def("get_self_transition", &Nonhomogeneous_markov::get_self_transition,
		return_value_policy<manage_new_object> (),
		python::args("index"),"return Function* of self transition")
	//.def("get_process", &Nonhomogeneous_markov::get_process,
	//	"return non parametric sequence process")

/*

    Nonhomogeneous_markov();
    Nonhomogeneous_markov(int inb_state , int *ident);
    Nonhomogeneous_markov(const Chain *pchain , const Function **pself_transition , int length);
    Nonhomogeneous_markov(const Nonhomogeneous_markov &markov , bool data_flag = true ,
                          bool characteristic_flag = true)
    :Chain(markov) { copy(markov , data_flag , characteristic_flag); }
    virtual ~Nonhomogeneous_markov();
    Nonhomogeneous_markov& operator=(const Nonhomogeneous_markov &markov);


    void characteristic_computation(int length , bool counting_flag);
    void characteristic_computation(const Nonhomogeneous_markov_data &seq , bool counting_flag ,
                                    bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq ,
                                  int index = I_DEFAULT) const;

    Nonhomogeneous_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                           bool counting_flag = true) const;
    Nonhomogeneous_markov_data* simulation(Format_error &error , int nb_sequence , int length ,
                                           bool counting_flag = true) const;
    Nonhomogeneous_markov_data* simulation(Format_error &error , int nb_sequence ,
                                           const Markovian_sequences &iseq ,
                                           bool counting_flag = true) const;

    Nonhomogeneous_markov_data* get_markov_data() const { return markov_data; }
    Nonparametric_sequence_process* get_process() const { return process; }
   */



;
}

