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



#ifndef HIDDEN_SEMI_MARKOV_H
#define HIDDEN_SEMI_MARKOV_H



/****************************************************************
 *
 *  Constantes :
 */


const double SEMI_MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // seuil pour stopper les iterations EM
const int EXPLORATION_NB_ITER = 10;    // nombre d'iterations de la phase d'exploration
const int STOCHASTIC_EXPLORATION_NB_ITER = 5;  // nombre d'iterations de la phase d'exploration
const int SEMI_MARKOV_NB_ITER = 500;   // nombre maximum d'iterations EM

const double MIN_SMOOTHED_PROBABILITY = 1.e-3;  // seuil sur les probabilites lissees

enum {
  SSTATE ,                             // etat
  IN_STATE ,                           // entree etat
  OUT_STATE                            // sortie etat
};



/****************************************************************
 *
 *  Definition des classes :
 */


class Hidden_semi_markov : public Semi_markov {  // semi-chaine de Markov cachee

    friend class Markovian_sequences;

    friend Hidden_semi_markov* hidden_semi_markov_ascii_read(Format_error &error , const char *path ,
                                                             int length,
                                                             bool counting_flag ,
                                                             double cumul_threshold ,
                                                             bool old_format );
    friend std::ostream& operator<<(std::ostream &os , const Hidden_semi_markov &hsmarkov)
    { return hsmarkov.ascii_write(os); }

private :

    Hidden_semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value)
    :Semi_markov(itype , inb_state , inb_output_process , nb_value) {}

    double forward_backward(const Markovian_sequences &seq , int index , std::ostream &os ,
                            int output , char format , double &max_marginal_entropy ,
                            double &entropy1) const;
    double forward_backward_sampling(const Markovian_sequences &seq , int index ,
                                     std::ostream &os , char format = 'a' ,
                                     int nb_state_sequence = NB_STATE_SEQUENCE) const;

    void log_computation();
    double viterbi(const Markovian_sequences &seq , double *posterior_probability = 0 ,
                   int index = I_DEFAULT) const;
    double generalized_viterbi(const Markovian_sequences &seq , int index , std::ostream &os ,
                               double seq_likelihood , char format , int inb_state_sequence) const;
    double viterbi_forward_backward(const Markovian_sequences &seq , int index ,
                                    std::ostream &os , int output , char format ,
                                    double seq_likelihood = D_INF) const;

public :

    Hidden_semi_markov() {}
    Hidden_semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                       int inb_output_process , Nonparametric_process **pobservation ,
                       int length , bool counting_flag)
    :Semi_markov(pchain , poccupancy , inb_output_process , pobservation , length ,
                 counting_flag) {}
    Hidden_semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                       int inb_output_process , Nonparametric_process **nonparametric_observation ,
                       Parametric_process **parametric_observation , int length , bool counting_flag)
    :Semi_markov(pchain , poccupancy , inb_output_process , nonparametric_observation ,
                 parametric_observation , length , counting_flag) {}
    Hidden_semi_markov(const Hidden_semi_markov &hsmarkov , bool data_flag = true ,
                       int param = I_DEFAULT)
    :Semi_markov(hsmarkov , data_flag , param) {}
    ~Hidden_semi_markov();

    Hidden_semi_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;

/*    RWDECLARE_COLLECTABLE(Hidden_semi_markov);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

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
                                                 bool characteristic_flag = false) const;

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
                                            const Markovian_sequences **seq , const char *path = 0) const;
};


Hidden_semi_markov* hidden_semi_markov_ascii_read(Format_error &error , const char *path ,
                                                  int length = DEFAULT_LENGTH ,
                                                  bool counting_flag = true ,
                                                  double cumul_threshold = OCCUPANCY_THRESHOLD ,
                                                  bool old_format = false);



#endif
