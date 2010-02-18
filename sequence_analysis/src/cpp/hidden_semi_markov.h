/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2010 CIRAD/INRIA Virtual Plants
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for V-Plants developers:
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


class HiddenSemiMarkov : public SemiMarkov {  // semi-chaine de Markov cachee

    friend class MarkovianSequences;

    friend HiddenSemiMarkov* hidden_semi_markov_ascii_read(StatError &error , const char *path ,
                                                           int length,
                                                           bool counting_flag ,
                                                           double cumul_threshold ,
                                                           bool old_format );
    friend std::ostream& operator<<(std::ostream &os , const HiddenSemiMarkov &hsmarkov)
    { return hsmarkov.ascii_write(os); }

private :

    HiddenSemiMarkov(char itype , int inb_state , int inb_output_process , int *nb_value)
    :SemiMarkov(itype , inb_state , inb_output_process , nb_value) {}

    int end_state() const;

    double forward_backward(const MarkovianSequences &seq , int index , std::ostream *os ,
                            MultiPlotSet *plot_set , int output , char format ,
                            double &max_marginal_entropy , double &entropy1) const;
    double forward_backward_sampling(const MarkovianSequences &seq , int index ,
                                     std::ostream &os , char format = 'a' ,
                                     int nb_state_sequence = NB_STATE_SEQUENCE) const;

    void log_computation();
    double viterbi(const MarkovianSequences &seq , double *posterior_probability = NULL ,
                   int index = I_DEFAULT) const;
    double generalized_viterbi(const MarkovianSequences &seq , int index , std::ostream &os ,
                               double seq_likelihood , char format , int inb_state_sequence) const;
    double viterbi_forward_backward(const MarkovianSequences &seq , int index ,
                                    std::ostream *os , MultiPlot *plot , int output ,
                                    char format , double seq_likelihood = D_INF) const;

public :

    HiddenSemiMarkov() {}
    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
                     int inb_output_process , NonparametricProcess **pobservation ,
                     int length , bool counting_flag)
    :SemiMarkov(pchain , poccupancy , inb_output_process , pobservation , length ,
                counting_flag) {}
    HiddenSemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
                     int inb_output_process , NonparametricProcess **nonparametric_observation ,
                     DiscreteParametricProcess **parametric_observation , int length , bool counting_flag)
    :SemiMarkov(pchain , poccupancy , inb_output_process , nonparametric_observation ,
                parametric_observation , length , counting_flag) {}
    HiddenSemiMarkov(const HiddenSemiMarkov &hsmarkov , bool data_flag = true ,
                     int param = I_DEFAULT)
    :SemiMarkov(hsmarkov , data_flag , param) {}
    ~HiddenSemiMarkov();

    HiddenSemiMarkov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;

    double likelihood_computation(const MarkovianSequences &seq , double *posterior_probability = NULL ,
                                  int index = I_DEFAULT) const;

    bool state_profile_write(StatError &error , std::ostream &os , const MarkovianSequences &iseq ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const char *path , const MarkovianSequences &iseq ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_ascii_write(StatError &error , std::ostream &os , int identifier ,
                                   int output = SSTATE , int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const char *path ,
                             int identifier = I_DEFAULT , int output = SSTATE ,
                             char format = 'a' , int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_plot_write(StatError &error , const char *prefix ,
                                  const MarkovianSequences &iseq , int identifier ,
                                  int output = SSTATE , const char *title = NULL) const;
    bool state_profile_plot_write(StatError &error , const char *prefix , int identifier ,
                                  int output = SSTATE , const char *title = NULL) const;

    MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                               const MarkovianSequences &iseq ,
                                               int identifier , int output = SSTATE) const;
    MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                               int identifier , int output = SSTATE) const;

    SemiMarkovData* state_sequence_computation(StatError &error , ostream &os ,
                                               const MarkovianSequences &seq ,
                                               bool characteristic_flag = true) const;

    SemiMarkovData* simulation(StatError &error , const FrequencyDistribution &hlength ,
                               bool counting_flag = true , bool divergence_flag = false) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence ,
                               int length , bool counting_flag = true) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence ,
                               const MarkovianSequences &iseq , bool counting_flag = true) const;

    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const HiddenSemiMarkov **ihsmarkov ,
                                           FrequencyDistribution **hlength ,
                                           const char *path = NULL) const;
    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const HiddenSemiMarkov **hsmarkov , int nb_sequence ,
                                           int length , const char *path = NULL) const;
    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const HiddenSemiMarkov **hsmarkov , int nb_sequence ,
                                           const MarkovianSequences **seq , const char *path = NULL) const;
};


HiddenSemiMarkov* hidden_semi_markov_ascii_read(StatError &error , const char *path ,
                                                int length = DEFAULT_LENGTH ,
                                                bool counting_flag = true ,
                                                double cumul_threshold = OCCUPANCY_THRESHOLD ,
                                                bool old_format = false);



#endif
