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
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/sequence_label.h"
#include "tool/config.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
//using namespace stat_tool;

class HiddenSemiMarkovWrap {

public:


};

// Boost declaration

void class_hidden_semi_markov() {

  enum_<sequence_analysis::wrap_util::UniqueInt<3, 102> >("SemiMarkovState")
    .value("SSTATE", SSTATE)
    .value("IN_STATE", IN_STATE)
    .value("OUT_STATE", OUT_STATE)
    .export_values();

	class_<Hidden_semi_markov, bases<Semi_markov > >
	("_Hidden_semi_markov", "Hidden_semi_markov")
    .def(init <Hidden_semi_markov, bool, int>())


/*
 *
    Hidden_semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                       int inb_output_process , Nonparametric_process **pobservation ,
                       int length , bool counting_flag)    :Semi_markov(pchain , poccupancy , inb_output_process , pobservation , length ,  counting_flag) {}
    Hidden_semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,    int inb_output_process , Nonparametric_process **nonparametric_observation ,
                       Parametric_process **parametric_observation , int length , bool counting_flag)    :Semi_markov(pchain , poccupancy , inb_output_process , nonparametric_observation ,                 parametric_observation , length , counting_flag) {}

    Hidden_semi_markov* thresholding(double min_probability = MIN_PROBABILITY) const;


    double likelihood_computation(const Markovian_sequences &seq , double *posterior_probability = 0 ,
                                  int index = I_DEFAULT) const;

    bool state_profile_write(Format_error &error , std::ostream &os , const Markovian_sequences &iseq ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(Format_error &error , const char *path , const Markovian_sequences &iseq ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_ascii_write(Format_error &error , std::ostream &os , int identifier ,
                                   int output = SSTATE , int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(Format_error &error , const char *path ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                  const Markovian_sequences &iseq , int identifier ,
                                  int output = SSTATE , const char *title = 0) const;
    bool state_profile_plot_write(Format_error &error , const char *prefix , int identifier ,
                                  int output = SSTATE , const char *title = 0) const;

    Semi_markov_data* state_sequence_computation(Format_error &error , ostream &os ,
                                                 const Markovian_sequences &seq ,
                                                 bool characteristic_flag = true) const;

    Semi_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                 bool counting_flag = true , bool divergence_flag = false) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence ,
                                 int length , bool counting_flag = true) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence ,
                                 const Markovian_sequences &iseq , bool counting_flag = true) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_semi_markov **ihsmarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_semi_markov **hsmarkov , int nb_sequence ,
                                            int length , const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_semi_markov **hsmarkov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;*/
;

}



