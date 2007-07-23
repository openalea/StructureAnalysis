/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#ifndef HIDDEN_VARIABLE_ORDER_MARKOV_H
#define HIDDEN_VARIABLE_ORDER_MARKOV_H



/****************************************************************
 *
 *  Constantes :
 */


const double VARIABLE_ORDER_MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // seuil pour stopper les iterations EM
const int VARIABLE_ORDER_MARKOV_NB_ITER = 100;  // nombre maximum d'iterations EM



/****************************************************************
 *
 *  Definition des classes :
 */


class Hidden_variable_order_markov : public Variable_order_markov {  // chaine de Markov
                                                                     // d'ordre variable cachee
    friend class Markovian_sequences;

    friend Hidden_variable_order_markov* hidden_variable_order_markov_ascii_read(Format_error &error ,
                                                                                 const char *path ,
                                                                                 int length,
                                                                                 double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const Hidden_variable_order_markov &hmarkov)
    { return hmarkov.ascii_write(os); }

private :

    double forward_backward(const Variable_order_markov_data &seq , int index ,
                            std::ostream &os , char format , double &max_marginal_entropy ,
                            double &entropy1) const;
    double forward_backward_sampling(const Variable_order_markov_data &seq , int index ,
                                     std::ostream &os , char format = 'a' ,
                                     int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_write(Format_error &error , std::ostream &os ,
                             const Variable_order_markov_data &iseq ,
                             int identifier = I_DEFAULT , char format = 'a' ,
                             int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                  const Variable_order_markov_data &iseq ,
                                  int identifier , const char *title = 0) const;

    void log_computation();
    double viterbi(const Variable_order_markov_data &seq , double *posterior_probability = 0 ,
                   int index = I_DEFAULT) const;
    double generalized_viterbi(const Variable_order_markov_data &seq , int index ,
                               std::ostream &os , double seq_likelihood , char format ,
                               int inb_state_sequence) const;
    double viterbi_forward_backward(const Variable_order_markov_data &seq , int index ,
                                    std::ostream &os , char format ,
                                    double seq_likelihood = D_INF) const;

public :

    Hidden_variable_order_markov() {}
    Hidden_variable_order_markov(const Variable_order_markov *pmarkov , int inb_output_process ,
                                 Nonparametric_process **nonparametric_observation ,
                                 Parametric_process **parametric_observation , int length)
    :Variable_order_markov(pmarkov , inb_output_process , nonparametric_observation ,
                           parametric_observation , length) {}
    Hidden_variable_order_markov(const Hidden_variable_order_markov &hmarkov , bool data_flag = true)
    :Variable_order_markov(hmarkov , data_flag) {}
    ~Hidden_variable_order_markov();

    Hidden_variable_order_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;

    double likelihood_computation(const Markovian_sequences &seq , double *posterior_probability = 0 ,
                                  int index = I_DEFAULT) const;

    bool state_profile_ascii_write(Format_error &error , std::ostream &os , int identifier ,
                                   int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(Format_error &error , const char *path ,
                             int identifier = I_DEFAULT , char format = 'a' ,
                             int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                  int identifier , const char *title = 0) const;

    Variable_order_markov_data* state_sequence_computation(Format_error &error ,
                                                           const Markovian_sequences &seq ,
                                                           bool characteristic_flag = false) const;

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
};


Hidden_variable_order_markov* hidden_variable_order_markov_ascii_read(Format_error &error ,
                                                                      const char *path ,
                                                                      int length = DEFAULT_LENGTH ,
                                                                      double cumul_threshold = OCCUPANCY_THRESHOLD);



#endif
