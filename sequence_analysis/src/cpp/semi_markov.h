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



#ifndef SEMI_MARKOV_H
#define SEMI_MARKOV_H



/****************************************************************
 *
 *  Constantes :
 */


const int LEAVE_LENGTH = 10000;         // longueur maximum pour le calcul de
                                        // la probabilite de quitter un etat

const double OCCUPANCY_LIKELIHOOD_DIFF = 1.e-5;  // seuil pour stopper les iterations EM
const int OCCUPANCY_NB_ITER = 10000;   // nombre maximum d'iterations EM
const int OCCUPANCY_COEFF = 10;        // coefficient arrondi estimateur pour les lois
                                       // d'occupation des etats

enum {
  MARKOVIAN ,
  SEMI_MARKOVIAN
};



/****************************************************************
 *
 *  Definition des classes :
 */


// class SemiMarkov : public StatInterface , public Chain {
class SemiMarkov : public StatInterface , protected Chain {  // semi-chaine de Markov

    friend class MarkovianSequences;
    friend class SemiMarkovIterator;
    friend class SemiMarkovData;

    friend SemiMarkov* semi_markov_ascii_read(StatError &error , const char *path ,
                                              int length, bool counting_flag ,
                                              double cumul_threshold );
    friend std::ostream& operator<<(std::ostream &os , const SemiMarkov &smarkov)
    { return smarkov.ascii_write(os , smarkov.semi_markov_data); }

protected :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    SemiMarkovData *semi_markov_data;  // pointeur sur un objet SemiMarkovData
    int *state_subtype;     //  MARKOVIAN/SEMI_MARKOVIAN
    Forward **forward;      // lois de l'intervalle de temps residuel
    int nb_output_process;  // nombre de processus d'observation
    NonparametricSequenceProcess **nonparametric_process;  // processus d'observation discrets non-parametriques
    DiscreteParametricProcess **discrete_parametric_process;  // processus d'observation discrets parametriques
    ContinuousParametricProcess **continuous_parametric_process;  // processus d'observation continus parametriques

    SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
               int inb_output_process , NonparametricProcess **pobservation ,
               int length , bool counting_flag);
    SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
               int inb_output_process , NonparametricProcess **nonparametric_observation ,
               DiscreteParametricProcess **discrete_parametric_observation ,
               ContinuousParametricProcess **continuous_parametric_observation ,
               int length , bool counting_flag);

    void copy(const SemiMarkov &smarkov , bool data_flag = true ,
              int param = I_DEFAULT);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const SemiMarkovData *seq ,
                              bool exhaustive = false , bool file_flag  = false ,
                              bool hidden = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const SemiMarkovData *seq ,
                                    bool hidden = false) const;
    bool plot_write(const char *prefix , const char *title ,
                    const SemiMarkovData *seq) const;
    MultiPlotSet* get_plotable(const SemiMarkovData *seq) const;

    int nb_parameter_computation(double min_probability = 0.) const;
    double penalty_computation(bool hidden , double min_probability = 0.) const;

    void initial_probability_computation();

    void index_state_distribution();
    double* memory_computation() const;
    void state_no_occurrence_probability(int state , double increment = LEAVE_INCREMENT);
    void state_first_occurrence_distribution(int state , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void state_leave_probability(int state , double increment = LEAVE_INCREMENT);
    void state_recurrence_time_distribution(int state , int min_nb_value = 1 ,
                                            double cumul_threshold = OCCUPANCY_THRESHOLD);
    void state_nb_pattern_mixture(int state , char pattern);

    void index_output_distribution(int variable);
    void output_no_occurrence_probability(int variable , int output ,
                                          double increment = LEAVE_INCREMENT);
    void output_first_occurrence_distribution(int variable , int output ,
                                              int min_nb_value = 1 ,
                                              double cumul_threshold = CUMUL_THRESHOLD);
    void output_leave_probability(const double *memory ,
                                  int variable , int output ,
                                  double increment = LEAVE_INCREMENT);
    void output_recurrence_time_distribution(const double *memory , int variable ,
                                             int output , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void output_sojourn_time_distribution(const double *memory , int variable ,
                                          int output , int min_nb_value = 1 ,
                                          double cumul_threshold = CUMUL_THRESHOLD);
    void output_nb_run_mixture(int variable , int output);
    void output_nb_occurrence_mixture(int variable , int output);

public :

    SemiMarkov();
    SemiMarkov(char itype , int inb_state , int inb_output_process , int *nb_value);
    SemiMarkov(const Chain *pchain , const NonparametricSequenceProcess *poccupancy ,
               const NonparametricProcess *pobservation , int length ,
               bool counting_flag);
    SemiMarkov(const SemiMarkov &smarkov , bool data_flag = true ,
               int param = I_DEFAULT)
    :Chain(smarkov) { copy(smarkov , data_flag , param); }
    void conditional_delete();
    ~SemiMarkov();
    SemiMarkov& operator=(const SemiMarkov &smarkov);

    DiscreteParametricModel* extract(StatError &error , int type ,
                                     int variable , int value) const;
    DiscreteParametricModel* extract(StatError &error , int state ,
                                     int frequency_distribution_type = FINAL_RUN) const;
    SemiMarkovData* extract_data(StatError &error) const;

    SemiMarkov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const SemiMarkovData &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    double likelihood_computation(const MarkovianSequences &seq , int index) const;
    double likelihood_computation(const SemiMarkovData &seq) const;

    SemiMarkovData* simulation(StatError &error , const FrequencyDistribution &hlength ,
                               bool counting_flag = true , bool divergence_flag = false) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence , int length ,
                               bool counting_flag = true) const;
    SemiMarkovData* simulation(StatError &error , int nb_sequence ,
                               const MarkovianSequences &iseq , bool counting_flag = true) const;

    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const SemiMarkov **ismarkov ,
                                           FrequencyDistribution **hlength ,
                                           const char *path = NULL) const;
    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const SemiMarkov **smarkov , int nb_sequence ,
                                           int length , const char *path = NULL) const;
    DistanceMatrix* divergence_computation(StatError &error , std::ostream &os , int nb_model ,
                                           const SemiMarkov **smarkov , int nb_sequence ,
                                           const MarkovianSequences **seq , const char *path = NULL) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    SemiMarkovData* get_semi_markov_data() const { return semi_markov_data; }
    int get_state_subtype(int state) const { return state_subtype[state]; }
    Forward** get_forward() const { return forward; }
    Forward* get_forward(int state) const { return forward[state]; }
    int get_nb_output_process() const { return nb_output_process; }
    NonparametricSequenceProcess* get_nonparametric_process(int variable)
    const { return nonparametric_process[variable]; }
    DiscreteParametricProcess** get_discrete_parametric_process() const
    { return discrete_parametric_process; }
    DiscreteParametricProcess* get_discrete_parametric_process(int variable) const
    { return discrete_parametric_process[variable]; }
    ContinuousParametricProcess** get_continuous_parametric_process() const
    { return continuous_parametric_process; }
    ContinuousParametricProcess* get_continuous_parametric_process(int variable) const
    { return continuous_parametric_process[variable]; }
};


SemiMarkov* semi_markov_ascii_read(StatError &error , const char *path ,
                                   int length = DEFAULT_LENGTH , bool counting_flag = true ,
                                   double cumul_threshold = OCCUPANCY_THRESHOLD);



class SemiMarkovIterator {  // iterateur semi-chaine de Markov

private :

    SemiMarkov *semi_markov;  // pointeur sur un objet SemiMarkov
    int state;              // etat
    int occupancy;          // temps d'occupation de l'etat
    int counter;            // compteur

    void copy(const SemiMarkovIterator &it);

public :

    SemiMarkovIterator(SemiMarkov *ismarkov);
    SemiMarkovIterator(const SemiMarkovIterator &iter)
    { copy(iter); }
    ~SemiMarkovIterator();
    SemiMarkovIterator& operator=(const SemiMarkovIterator &iter);

    bool simulation(int **int_seq , int length = 1 , bool initialization = false);
    int** simulation(int length = 1 , bool initialization = false);

    // acces membres de la classe

    SemiMarkov* get_semi_markov() const { return semi_markov; }
    int get_state() const { return state; }
    int get_occupancy() const { return occupancy; }
    int get_counter() const { return counter; }
    int get_nb_variable() const { return (semi_markov ? semi_markov->nb_output_process + 1 : 0); }
};



class SemiMarkovData : public MarkovianSequences {  // structure de donnees correspondant
                                                    // a une semi-chaine de Markov

    friend class MarkovianSequences;
    friend class SemiMarkov;
    friend class HiddenSemiMarkov;

    friend std::ostream& operator<<(std::ostream &os , const SemiMarkovData &seq)
    { return seq.ascii_write(os , false); }

private :

    SemiMarkov *semi_markov;  // pointeur sur un objet SemiMarkov
    ChainData *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences
    double hidden_likelihood;  // vraisemblance de toutes les sequences possibles
    double sample_entropy;  // entropie des sequences d'etats
    double *posterior_probability;  // probabilite a posteriori de la sequence d'etats la plus probable
    double *entropy;        // entropie des sequences d'etats

    void copy(const SemiMarkovData &seq , bool model_flag = true);

public :

    SemiMarkovData();
    SemiMarkovData(const FrequencyDistribution &ihlength , int inb_variable ,
                   int *itype , bool init_flag = false);
    SemiMarkovData(const MarkovianSequences &seq);
    SemiMarkovData(const MarkovianSequences &seq , char transform , bool initial_run_flag);
    SemiMarkovData(const SemiMarkovData &seq , bool model_flag = true , char transform = 'c')
    :MarkovianSequences(seq , transform) { copy(seq , model_flag); }
    ~SemiMarkovData();
    SemiMarkovData& operator=(const SemiMarkovData &seq);

    DiscreteDistributionData* extract(StatError &error , int type ,
                                      int variable , int value) const;
    SemiMarkovData* remove_index_parameter(StatError &error) const;
    MarkovianSequences* build_auxiliary_variable(StatError &error) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(StatError &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(StatError &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(StatError &error , const char *path) const;
    bool plot_write(StatError &error , const char *prefix ,
                    const char *title = NULL) const;
    MultiPlotSet* get_plotable() const;

    void build_transition_count(const SemiMarkov *smarkov = NULL);

    // acces membres de la classe

    SemiMarkov* get_semi_markov() const { return semi_markov; }
    ChainData* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
    double get_hidden_likelihood() const { return hidden_likelihood; }
    double get_sample_entropy() const { return sample_entropy; }
    double get_posterior_probability(int index) const { return posterior_probability[index]; }
    double get_entropy(int index) const { return entropy[index]; }
};



#endif
