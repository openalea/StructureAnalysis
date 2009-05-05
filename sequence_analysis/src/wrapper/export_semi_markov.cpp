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

using namespace boost::python;
using namespace boost;
//using namespace stat_tool;

class SemiMarkovWrap {

public:


};

// Boost declaration

void class_semi_markov() {

	class_<Semi_markov, bases<STAT_interface> >
	("_Semi_markov", "Semi_markov")
;



/*
 *
  Semi_markov();
    Semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value);
    Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                const Nonparametric_process *pobservation , int length ,
                bool counting_flag);
    Semi_markov(const Semi_markov &smarkov , bool data_flag = true ,  int param = I_DEFAULT)  :Chain(smarkov) { copy(smarkov , data_flag , param); }

    Parametric_model* extract(Format_error &error , int type ,
                              int variable , int value) const;
    Parametric_model* extract(Format_error &error , int state ,
                              int histogram_type = FINAL_RUN) const;
    Semi_markov_data* extract_data(Format_error &error) const;

    Semi_markov* thresholding(double min_probability = MIN_PROBABILITY) const;



    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const Semi_markov_data &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq , int index) const;
    double likelihood_computation(const Semi_markov_data &seq) const;

    Semi_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                 bool counting_flag = true , bool divergence_flag = false) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence , int length ,
                                 bool counting_flag = true) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence ,
                                 const Markovian_sequences &iseq , bool counting_flag = true) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **ismarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **smarkov , int nb_sequence ,
                                            int length , const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **smarkov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;



int get_nb_iterator() const { return nb_iterator; }
Semi_markov_data* get_semi_markov_data() const { return semi_markov_data; }
int get_state_subtype(int state) const { return state_subtype[state]; }
Forward** get_forward() const { return forward; }
Forward* get_forward(int state) const { return forward[state]; }
int get_nb_output_process() const { return nb_output_process; }
Nonparametric_sequence_process* get_nonparametric_process(int variable)
const { return nonparametric_process[variable]; }
Parametric_process** get_parametric_process() const { return parametric_process; }
Parametric_process* get_parametric_process(int variable)
const { return parametric_process[variable]; }

*/

   ;
}


void class_semi_markov_data() {

	class_<Semi_markov_data, bases<Markovian_sequences > >
	("_Semi_markov_data", "Semi_markov_data")
	.def(init<Markovian_sequences>())
	.def(init<Markovian_sequences, char, bool>())
	.def(init<Semi_markov_data, bool, char>())

;


/*   Semi_markov_data();
Semi_markov_data(const Histogram &ihlength , int inb_variable , bool init_flag = false);

Distribution_data* extract(Format_error &error , int type ,
int variable , int value) const;
Semi_markov_data* remove_index_parameter(Format_error &error) const;


    void build_transition_count(const Semi_markov *smarkov = 0);
  
 Semi_markov* get_semi_markov() const { return semi_markov; }
 Chain_data* get_chain_data() const { return chain_data; }
 double get_likelihood() const { return likelihood; }
 double get_hidden_likelihood() const { return hidden_likelihood; }
 double get_posterior_probability(int index) const { return posterior_probability[index]; }
*/


}



