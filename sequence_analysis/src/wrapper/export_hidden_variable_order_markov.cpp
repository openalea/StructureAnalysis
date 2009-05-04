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
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/stat_label.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
#include "sequence_analysis/sequence_label.h"

#include "tool/config.h"

#include <boost/python.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/list.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/python/make_constructor.hpp>

using namespace boost::python;
using namespace boost;
using namespace stat_tool;








void class_hidden_variable_order_markov() {

        class_<Hidden_variable_order_markov, bases<Variable_order_markov > >
        ("_Hidden_variable_order_markov", "Hidden_variable_order_markov")
        //,def(init <const Hidden_variable_order_markov &, optional<bool> >())
;

      /*
        Hidden_variable_order_markov() {}
          Hidden_variable_order_markov(const Variable_order_markov *pmarkov , int inb_output_process ,
                                       Nonparametric_process **nonparametric_observation ,
                                       Parametric_process **parametric_observation , int length)


          Hidden_variable_order_markov(const Hidden_variable_order_markov &hmarkov , bool data_flag = true)


          Hidden_variable_order_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

          std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
          bool ascii_write(Format_error &error , const char *path ,
                           bool exhaustive = false) const;
          bool spreadsheet_write(Format_error &error , const char *path) const;

          double likelihood_computation(const Markovian_sequences &seq , double *posterior_probability = 0 ,
                                        int index = I_DEFAULT) const;

          bool state_profile_write(Format_error &error , std::ostream &os , const Markovian_sequences &iseq ,
                                   int identifier = I_DEFAULT , char format = 'a' ,
                                   int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(Format_error &error , const char *path , const Markovian_sequences &iseq ,
                                   int identifier = I_DEFAULT , char format = 'a' ,
                                   int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_ascii_write(Format_error &error , std::ostream &os , int identifier ,
                                         int state_sequence = GENERALIZED_VITERBI ,
                                         int nb_state_sequence = NB_STATE_SEQUENCE) const;
          bool state_profile_write(Format_error &error , const char *path ,
                                   int identifier = I_DEFAULT , char format = 'a' ,
                                   int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;

          bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                        const Markovian_sequences &iseq ,
                                        int identifier , const char *title = 0) const;
          bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                        int identifier , const char *title = 0) const;

          Variable_order_markov_data* state_sequence_computation(Format_error &error ,
                                                                 const Markovian_sequences &seq ,
                                                                 bool characteristic_flag = true) const;

          Variable_order_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                                 bool counting_flag = true , bool divergence_flag = false) const;
          Variable_order_markov_data* simulation(Format_error &error , int nb_sequence ,
                                                 int length , bool counting_flag = true) const;
          Variable_order_markov_data* simulation(Format_error &error , int nb_sequence ,
                                                 const Markovian_sequences &iseq ,
                                                 bool counting_flag = true) const;

          Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                                  const Hidden_variable_order_markov **ihmarkov ,
                                                  Histogram **hlength , const char *path = 0) const;
          Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                                  const Hidden_variable_order_markov **hmarkov , int nb_sequence ,
                                                  int length , const char *path = 0) const;
          Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                                  const Hidden_variable_order_markov **hmarkov , int nb_sequence ,
                                                  const Markovian_sequences **seq , const char *path = 0) const;

*/
}


