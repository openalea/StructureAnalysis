/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2015 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: hidden_variable_order_markov.h 18053 2015-04-23 09:44:23Z guedon $
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



#ifndef HIDDEN_VARIABLE_ORDER_MARKOV_H
#define HIDDEN_VARIABLE_ORDER_MARKOV_H



namespace sequence_analysis {



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


  class HiddenVariableOrderMarkov : public VariableOrderMarkov {  // chaine de Markov
                                                                  // d'ordre variable cachee
    friend class MarkovianSequences;

    friend HiddenVariableOrderMarkov* hidden_variable_order_markov_ascii_read(stat_tool::StatError &error ,
                                                                              const char *path ,
                                                                              int length,
                                                                              double cumul_threshold);
    friend std::ostream& operator<<(std::ostream &os , const HiddenVariableOrderMarkov &hmarkov)
    { return hmarkov.ascii_write(os); }

  private :

    void forward_backward(VariableOrderMarkovData &seq) const;
    double forward_backward(const MarkovianSequences &seq , int index ,
                            std::ostream *os , stat_tool::MultiPlotSet *plot_set , char format ,
                            double &max_marginal_entropy , double &entropy1) const;
    double forward_backward_sampling(const MarkovianSequences &seq , int index ,
                                     std::ostream &os , char format = 'a' ,
                                     int nb_state_sequence = NB_STATE_SEQUENCE) const;

    void log_computation();
    double viterbi(const MarkovianSequences &seq , double *posterior_probability ,
                   int index = I_DEFAULT) const;
    void viterbi(VariableOrderMarkovData &seq) const;
    double generalized_viterbi(const MarkovianSequences &seq , int index ,
                               std::ostream &os , double seq_likelihood , char format ,
                               int inb_state_sequence) const;
    double viterbi_forward_backward(const MarkovianSequences &seq , int index ,
                                    std::ostream *os , stat_tool::MultiPlot *plot , char format ,
                                    double seq_likelihood = D_INF) const;

  public :

    HiddenVariableOrderMarkov() {}
    HiddenVariableOrderMarkov(const VariableOrderMarkovChain *pmarkov , int inb_output_process ,
                              stat_tool::CategoricalProcess **categorical_observation ,
                              stat_tool::DiscreteParametricProcess **discrete_parametric_observation ,
                              stat_tool::ContinuousParametricProcess **continuous_parametric_observation ,
                              int length)
    :VariableOrderMarkov(pmarkov , inb_output_process , categorical_observation ,
                         discrete_parametric_observation ,
                         continuous_parametric_observation , length) {}
    HiddenVariableOrderMarkov(const HiddenVariableOrderMarkov &hmarkov , bool data_flag = true)
    :VariableOrderMarkov(hmarkov , data_flag) {}
    ~HiddenVariableOrderMarkov();

    HiddenVariableOrderMarkov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const char *path) const;

    double likelihood_computation(const MarkovianSequences &seq , double *posterior_probability = NULL ,
                                  int index = I_DEFAULT) const;

    bool state_profile_write(stat_tool::StatError &error , std::ostream &os ,
                             const MarkovianSequences &iseq ,
                             int identifier = I_DEFAULT , char format = 'a' ,
                             int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(stat_tool::StatError &error , const char *path ,
                             const MarkovianSequences &iseq ,
                             int identifier = I_DEFAULT , char format = 'a' ,
                             int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_ascii_write(stat_tool::StatError &error , std::ostream &os , int identifier ,
                                   int state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(stat_tool::StatError &error , const char *path ,
                             int identifier = I_DEFAULT , char format = 'a' ,
                             int state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_plot_write(stat_tool::StatError &error , const char *prefix ,
                                  const MarkovianSequences &iseq ,
                                  int identifier , const char *title = NULL) const;
    bool state_profile_plot_write(stat_tool::StatError &error , const char *prefix ,
                                  int identifier , const char *title = NULL) const;

    stat_tool::MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                                          const MarkovianSequences &iseq ,
                                                          int identifier) const;
    stat_tool::MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                                          int identifier) const;

    VariableOrderMarkovData* state_sequence_computation(stat_tool::StatError &error ,
                                                        const MarkovianSequences &seq ,
                                                        bool characteristic_flag = true) const;

    VariableOrderMarkovData* simulation(stat_tool::StatError &error ,
                                        const stat_tool::FrequencyDistribution &hlength ,
                                        bool counting_flag = true , bool divergence_flag = false) const;
    VariableOrderMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                        int length , bool counting_flag = true) const;
    VariableOrderMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                        const MarkovianSequences &iseq ,
                                        bool counting_flag = true) const;

    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const HiddenVariableOrderMarkov **ihmarkov ,
                                                      FrequencyDistribution **hlength ,
                                                      const char *path = NULL) const;
    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const HiddenVariableOrderMarkov **hmarkov , int nb_sequence ,
                                                      int length , const char *path = NULL) const;
    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const HiddenVariableOrderMarkov **hmarkov , int nb_sequence ,
                                                      const MarkovianSequences **seq , const char *path = NULL) const;
  };


  HiddenVariableOrderMarkov* hidden_variable_order_markov_ascii_read(stat_tool::StatError &error ,
                                                                     const char *path ,
                                                                     int length = DEFAULT_LENGTH ,
                                                                     double cumul_threshold = OCCUPANCY_THRESHOLD);


};  // namespace sequence_analysis



#endif
