/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       V-Plants: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2017 CIRAD/INRA/Inria Virtual Plants
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
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


#include "semi_markov.h"


namespace sequence_analysis {



/****************************************************************
 *
 *  Constants
 */


  const double SEMI_MARKOV_LIKELIHOOD_DIFF = 1.e-6;  // threshold for stopping the EM iterations
  const int EXPLORATION_NB_ITER = 10;    // number of iterations, exploration phase
  const int STOCHASTIC_EXPLORATION_NB_ITER = 5;  // number of iterations, exploration phase (MCEM algorithm)
  const int SEMI_MARKOV_NB_ITER = 500;   // maximum number of EM iterations

  const double MIN_SMOOTHED_PROBABILITY = 1.e-3;  // threshold on the smoothed probabilities

  enum state_profile {
    SSTATE ,                             // state
    IN_STATE ,                           // state entering
    OUT_STATE                            // state exit
  };



/****************************************************************
 *
 *  Class definition
 */


  /// \brief Hidden semi-Markov chain

  class HiddenSemiMarkov : public SemiMarkov {

    friend class MarkovianSequences;

    friend std::ostream& operator<<(std::ostream &os , const HiddenSemiMarkov &hsmarkov)
    { return hsmarkov.ascii_write(os); }

  private :

    HiddenSemiMarkov(stat_tool::process_type itype , int inb_state , int inb_output_process , int *nb_value)
    :SemiMarkov(itype , inb_state , inb_output_process , nb_value) {}

    int end_state() const;

    void forward_backward(SemiMarkovData &seq) const;
    double forward_backward(MarkovianSequences &seq , int index , std::ostream *os ,
                            stat_tool::MultiPlotSet *plot_set ,
                            state_profile output , stat_tool::output_format format ,
                            double &max_marginal_entropy , double &entropy1) const;
    double forward_backward_sampling(const MarkovianSequences &seq , int index , std::ostream &os ,
                                     stat_tool::output_format format = stat_tool::ASCII ,
                                     int nb_state_sequence = NB_STATE_SEQUENCE) const;

    void log_computation();
    double viterbi(const MarkovianSequences &seq , double *posterior_probability ,
                   int index = stat_tool::I_DEFAULT) const;
    void viterbi(SemiMarkovData &seq) const;
    double generalized_viterbi(const MarkovianSequences &seq , int index , std::ostream &os ,
                               double seq_likelihood , stat_tool::output_format format ,
                               int inb_state_sequence) const;
    double viterbi_forward_backward(const MarkovianSequences &seq , int index ,
                                    std::ostream *os , stat_tool::MultiPlot *plot ,
                                    state_profile output , stat_tool::output_format format ,
                                    double seq_likelihood = stat_tool::D_INF) const;

    bool state_profile_write(StatError &error , std::ostream &os , const MarkovianSequences &iseq ,
                             int identifier , state_profile output = SSTATE ,
                             stat_tool::output_format format = stat_tool::ASCII ,
                             latent_structure_algorithm state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;

  public :

    HiddenSemiMarkov() {}
    HiddenSemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                     int inb_output_process , stat_tool::CategoricalProcess **pobservation ,
                     int length , bool counting_flag)
    :SemiMarkov(pchain , poccupancy , inb_output_process , pobservation , length ,
                counting_flag) {}
    HiddenSemiMarkov(const Chain *pchain , const CategoricalSequenceProcess *poccupancy ,
                     int inb_output_process , stat_tool::CategoricalProcess **categorical_observation ,
                     stat_tool::DiscreteParametricProcess **discrete_parametric_observation ,
                     stat_tool::ContinuousParametricProcess **continuous_parametric_observation ,
                     int length , bool counting_flag)
    :SemiMarkov(pchain , poccupancy , inb_output_process , categorical_observation ,
                discrete_parametric_observation ,
                continuous_parametric_observation , length , counting_flag) {}
    HiddenSemiMarkov(const HiddenSemiMarkov &hsmarkov , bool data_flag = true ,
                     int param = stat_tool::I_DEFAULT)
    :SemiMarkov(hsmarkov , data_flag , param) {}
    ~HiddenSemiMarkov();

    HiddenSemiMarkov* thresholding(double min_probability = MIN_PROBABILITY) const;

    static HiddenSemiMarkov* ascii_read(StatError &error , const std::string path ,
                                        int length = DEFAULT_LENGTH ,
                                        bool counting_flag = true ,
                                        double cumul_threshold = OCCUPANCY_THRESHOLD ,
                                        bool old_format = false);

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const std::string path , bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const std::string path) const;

    double likelihood_computation(const MarkovianSequences &seq , double *posterior_probability = NULL ,
                                  int index = stat_tool::I_DEFAULT) const;

    bool state_profile_ascii_write(StatError &error , const MarkovianSequences &iseq ,
                                   int identifier , state_profile output = SSTATE ,
                                   latent_structure_algorithm state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const std::string path , const MarkovianSequences &iseq ,
                             int identifier , state_profile output = SSTATE ,
                             stat_tool::output_format format = stat_tool::ASCII ,
                             latent_structure_algorithm state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_ascii_write(StatError &error , int identifier , state_profile output = SSTATE ,
                                   latent_structure_algorithm state_sequence = GENERALIZED_VITERBI ,
                                   int nb_state_sequence = NB_STATE_SEQUENCE) const;
    bool state_profile_write(StatError &error , const std::string path ,
                             int identifier , state_profile output = SSTATE ,
                             stat_tool::output_format format = stat_tool::ASCII ,
                             latent_structure_algorithm state_sequence = GENERALIZED_VITERBI ,
                             int nb_state_sequence = NB_STATE_SEQUENCE) const;

    bool state_profile_plot_write(StatError &error , const char *prefix ,
                                  const MarkovianSequences &iseq , int identifier ,
                                  state_profile output = SSTATE , const char *title = NULL) const;
    bool state_profile_plot_write(StatError &error , const char *prefix , int identifier ,
                                  state_profile output = SSTATE , const char *title = NULL) const;

    stat_tool::MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                                          const MarkovianSequences &iseq ,
                                                          int identifier , state_profile output = SSTATE) const;
    stat_tool::MultiPlotSet* state_profile_plotable_write(StatError &error ,
                                                          int identifier , state_profile output = SSTATE) const;

    SemiMarkovData* state_sequence_computation(StatError &error , bool display ,
                                               const MarkovianSequences &seq ,
                                               bool characteristic_flag = true) const;

    SemiMarkovData* simulation(StatError &error , const FrequencyDistribution &hlength ,
                               bool counting_flag = true , bool divergence_flag = false) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence ,
                               int length , bool counting_flag = true) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence ,
                               const MarkovianSequences &iseq , bool counting_flag = true) const;

    stat_tool::DistanceMatrix* divergence_computation(StatError &error , bool display , int nb_model ,
                                                      const HiddenSemiMarkov **ihsmarkov ,
                                                      stat_tool::FrequencyDistribution **hlength ,
                                                      const std::string path = "") const;
    stat_tool::DistanceMatrix* divergence_computation(StatError &error , bool display , int nb_model ,
                                                      const HiddenSemiMarkov **hsmarkov , int nb_sequence ,
                                                      int length , const std::string path = "") const;
    stat_tool::DistanceMatrix* divergence_computation(StatError &error , bool display , int nb_model ,
                                                      const HiddenSemiMarkov **hsmarkov , int nb_sequence ,
                                                      const MarkovianSequences **seq , const std::string path = "") const;
  };


};  // namespace sequence_analysis



#endif
