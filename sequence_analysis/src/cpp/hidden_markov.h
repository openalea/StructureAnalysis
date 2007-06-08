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



#ifndef HIDDEN_MARKOV_H
#define HIDDEN_MARKOV_H



/****************************************************************
 *
 *  Constantes :
 */


const double MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // seuil pour stopper les iterations EM
const int MARKOV_NB_ITER = 100;        // nombre maximum d'iterations EM



/****************************************************************
 *
 *  Definition des classes :
 */


class Hidden_markov : public Markov {  // chaine de Markov cachee

    friend class Markovian_sequences;

    friend Hidden_markov* hidden_markov_ascii_read(Format_error &error , const char *path ,
                                                   int length,
                                                   bool old_format);
    friend std::ostream& operator<<(std::ostream &os , const Hidden_markov &hmarkov)
    { return hmarkov.ascii_write(os); }

private :

    double forward_backward(const Markov_data &seq , int index , std::ostream &os ,
                            char format) const;
    bool state_profile_write(Format_error &error , std::ostream &os , const Markov_data &iseq ,
                             int identifier = I_DEFAULT , char format = 'a') const;
    bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                  const Markov_data &iseq , int identifier ,
                                  const char *title = 0) const;

    void init(bool left_right , double self_transition);
    void init(double self_transition);
    void log_computation();
    double viterbi(const Markov_data &seq , int index = I_DEFAULT) const;

public :

    Hidden_markov() {}
    Hidden_markov(int inb_state , int iorder , int inb_output_process , int *nb_value)
    :Markov(inb_state , iorder , inb_output_process , nb_value) {}
    Hidden_markov(const Chain *pchain , int iorder , int inb_output_process ,
                  Nonparametric_process **pobservation , int length)
    :Markov(pchain , iorder , inb_output_process , pobservation , length) {}
    Hidden_markov(const Hidden_markov &hmarkov , bool data_flag = true ,
                  bool characteristic_flag = true)
    :Markov(hmarkov , data_flag , characteristic_flag) {}
    Hidden_markov(const Hidden_markov &hmarkov , double self_transition);
    Hidden_markov(const Hidden_markov &hmarkov , int state)
    :Markov(hmarkov , state) {}
    ~Hidden_markov();

    Hidden_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;

/*    RWDECLARE_COLLECTABLE(Hidden_markov);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    double likelihood_computation(const Markovian_sequences &seq , int index = I_DEFAULT) const;

    bool state_profile_ascii_write(Format_error &error , std::ostream &os , int identifier) const;
    bool state_profile_write(Format_error &error , const char *path ,
                             int identifier = I_DEFAULT , char format = 'a') const;
    bool state_profile_plot_write(Format_error &error , const char *prefix ,
                                  int identifier , const char *title = 0) const;

    Markov_data* state_sequence_computation(Format_error &error , const Markovian_sequences &seq ,
                                            bool characteristic_flag = false) const;

    Markov_data* simulation(Format_error &error , const Histogram &hlength ,
                            bool divergence_flag = false) const;
    Markov_data* simulation(Format_error &error , int nb_sequence , int length) const;
    Markov_data* simulation(Format_error &error , int nb_sequence ,
                            const Markovian_sequences &iseq) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_markov **ihmarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_markov **hmarkov , int nb_sequence ,
                                            int length , const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Hidden_markov **hmarkov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;
};


Hidden_markov* hidden_markov_ascii_read(Format_error &error , const char *path ,
                                        int length = DEFAULT_LENGTH ,
                                        bool old_format = false);



#endif
