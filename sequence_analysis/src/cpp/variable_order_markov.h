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
 *       $Id: variable_order_markov.h 18081 2015-04-23 11:00:24Z guedon $
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



#ifndef VARIABLE_ORDER_MARKOV_H
#define VARIABLE_ORDER_MARKOV_H



namespace Stat_trees {
  class MarkovOutTree;
  MarkovOutTree* markov_out_tree_parsing(StatError& error,
                                         std::ifstream &in_file,
                                         int &line);
};


namespace sequence_analysis {



/****************************************************************
 *
 *  Constantes :
 */


  const int MAX_LAG = 100;               // decalage maximum pour le calcul des coefficients
                                         // d'autocorrelation
  const int MEMORY_MIN_COUNT = 10;       // effectif minimum pour comparer
                                         // une memoire a ses fils
  const double LAPLACE_COEFF = 1.;       // coefficient pour l'estimateur de Laplace

  enum {
    NON_TERMINAL ,
    TERMINAL ,
    COMPLETION
  };



/****************************************************************
 *
 *  Definition des classes :
 */


  class VariableOrderMarkovChain : public stat_tool::Chain {  // chaine de Markov d'ordre variable

  public :

    int *memory_type;       // types des memoires (NON_TERMINAL/TERMINAL/COMPLETION)
    int *order;             // ordres des memoires
    int max_order;          // ordre maximum des memoires
    int **state;            // etats composant les memoires
    int *parent;            // memoires parent
    int **child;            // memoires fils
    int **next;             // memoires suivantes
    int *nb_memory;         // nombre de memoires precedentes
    int **previous;         // memoires precedentes
    CategoricalSequenceProcess *state_process;  // processus d'etat

    void memory_tree_completion(const VariableOrderMarkovChain &markov);
    void build(const VariableOrderMarkovChain &markov);
    void copy(const VariableOrderMarkovChain &markov);
    void remove();

    VariableOrderMarkovChain();
    VariableOrderMarkovChain(char itype , int inb_state , int inb_row);
    VariableOrderMarkovChain(char itype , int inb_state , int inb_row , int imax_order);
    VariableOrderMarkovChain(char itype , int inb_state , int iorder , bool init_flag);
    VariableOrderMarkovChain(const VariableOrderMarkovChain &markov)
    :Chain(markov) { copy(markov); }
    ~VariableOrderMarkovChain();
    VariableOrderMarkovChain& operator=(const VariableOrderMarkovChain &markov);

    void find_parent_memory(int index);
    void build_memory_transition();
    void build_previous_memory();
    bool check_free_suffix() const;
    bool** logic_transition_computation() const;
    void component_computation();

    void build_non_terminal();

    void thresholding(double min_probability);

    void max_order_computation();
    int nb_parameter_computation(double min_probability = 0.) const;
    int nb_transient_parameter_computation(double min_probability = 0.) const;

    std::ostream& ascii_memory_tree_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& ascii_transition_tree_print(std::ostream &os , bool file_flag = false) const;

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;

    void non_terminal_transition_probability_computation();
    void initial_probability_computation();

    void index_state_distribution();
    double* memory_computation() const;
    void state_no_occurrence_probability(int istate , double increment = LEAVE_INCREMENT);
    void state_first_occurrence_distribution(int istate , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void state_leave_probability(const double *imemory , int istate ,
                                 double increment = LEAVE_INCREMENT);
    void state_recurrence_time_distribution(const double *imemory , int istate ,
                                            int min_nb_value = 1 ,
                                            double cumul_threshold = CUMUL_THRESHOLD);
    void state_sojourn_time_distribution(const double *imemory , int istate ,
                                         int min_nb_value = 1 ,
                                         double cumul_threshold = CUMUL_THRESHOLD);
    void state_nb_pattern_mixture(int istate , char pattern);

    Correlation* state_autocorrelation_computation(stat_tool::StatError &error ,
                                                   int istate , int max_lag ,
                                                   const MarkovianSequences *seq) const;
  };


  VariableOrderMarkovChain* variable_order_markov_parsing(stat_tool::StatError &error ,
                                                          std::ifstream &in_file ,
                                                          int &line , char type);



  class VariableOrderMarkovChainData;
  class VariableOrderMarkovData;


  class VariableOrderMarkov : public stat_tool::StatInterface , protected VariableOrderMarkovChain {  // chaine de Markov
                                                                                                      // d'ordre variable
    friend class MarkovianSequences;
    friend class VariableOrderMarkovIterator;
    friend class VariableOrderMarkovChainData;
    friend class VariableOrderMarkovData;
    friend class Stat_trees::MarkovOutTree;  // revoir avec J.B.

    friend VariableOrderMarkov* variable_order_markov_ascii_read(stat_tool::StatError &error ,
                                                                 const char *path ,int length);
    friend Stat_trees::MarkovOutTree* Stat_trees::markov_out_tree_parsing(stat_tool::StatError& error,  // revoir avec J.B.
                                                                          std::ifstream &in_file, int &line);
                                                             
    friend std::ostream& operator<<(std::ostream &os , const VariableOrderMarkov &markov)
    { return markov.ascii_write(os); }

  protected :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    VariableOrderMarkovData *markov_data;  // pointeur sur un objet VariableOrderMarkovData
    int nb_output_process;  // nombre de processus d'observation
    CategoricalSequenceProcess **categorical_process;  // processus d'observation categoriels
    stat_tool::DiscreteParametricProcess **discrete_parametric_process;  // processus d'observation discrets parametriques
    stat_tool::ContinuousParametricProcess **continuous_parametric_process;  // processus d'observation continus parametriques

    VariableOrderMarkov(const VariableOrderMarkovChain *pmarkov , int inb_output_process ,
                        stat_tool::CategoricalProcess **categorical_observation ,
                        stat_tool::DiscreteParametricProcess **discrete_parametric_observation ,
                        stat_tool::ContinuousParametricProcess **continuous_parametric_observation ,
                        int length);

    void copy(const VariableOrderMarkov &markov , bool data_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const VariableOrderMarkovData *seq ,
                              bool exhaustive = false , bool file_flag = false ,
                              bool hidden = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const VariableOrderMarkovData *seq ,
                                    bool hidden = false) const;
    bool plot_write(const char *prefix , const char *title ,
                    const VariableOrderMarkovData *seq) const;
    stat_tool::MultiPlotSet* get_plotable(const VariableOrderMarkovData *seq) const;

    int nb_parameter_computation(double min_probability = 0.) const;
    double penalty_computation(bool hidden , double min_probability = 0.) const;

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

    Correlation* output_autocorrelation_computation(stat_tool::StatError &error , int variable ,
                                                    int output , int max_lag ,
                                                    const VariableOrderMarkovData *seq) const;

    double likelihood_computation(const VariableOrderMarkovChainData &chain_data) const;

    double likelihood_correction(const VariableOrderMarkovData &seq) const;

    std::ostream& transition_count_ascii_write(std::ostream &os , bool begin) const;

  public :

    VariableOrderMarkov();
    VariableOrderMarkov(char itype , int inb_state , int inb_row);
    VariableOrderMarkov(char itype , int inb_state , int inb_row , int imax_order);
    VariableOrderMarkov(char itype , int inb_state , int iorder , bool init_flag ,
                        int inb_output_process = 0 , int nb_value = 0);
    VariableOrderMarkov(const VariableOrderMarkov &markov ,
                        int inb_output_process , int nb_value);
/*    VariableOrderMarkov(const VariableOrderMarkov &markov ,
                        int inb_output_process , int *nb_value); */
    VariableOrderMarkov(const VariableOrderMarkovChain *pmarkov ,
                        const stat_tool::CategoricalProcess *pobservation , int length);
    VariableOrderMarkov(const VariableOrderMarkov &markov , bool data_flag = true)
    :VariableOrderMarkovChain(markov) { copy(markov , data_flag); }
    void conditional_delete();
    ~VariableOrderMarkov();
    VariableOrderMarkov& operator=(const VariableOrderMarkov &markov);

    DiscreteParametricModel* extract(stat_tool::StatError &error , int type ,
                                     int variable , int value) const;
    VariableOrderMarkovData* extract_data(stat_tool::StatError &error) const;

    VariableOrderMarkov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const char *path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const VariableOrderMarkovData &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    Correlation* state_autocorrelation_computation(stat_tool::StatError &error , int istate ,
                                                   int max_lag = MAX_LAG) const;
    Correlation* output_autocorrelation_computation(stat_tool::StatError &error , int variable ,
                                                    int output , int max_lag = MAX_LAG) const;

    double likelihood_computation(const MarkovianSequences &seq , int index) const;
    double likelihood_computation(const VariableOrderMarkovData &seq) const;

    VariableOrderMarkovData* simulation(stat_tool::StatError &error , const FrequencyDistribution &hlength ,
                                        bool counting_flag = true , bool divergence_flag = false) const;
    VariableOrderMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                        int length , bool counting_flag = true) const;
    VariableOrderMarkovData* simulation(stat_tool::StatError &error , int nb_sequence ,
                                        const MarkovianSequences &iseq ,
                                        bool counting_flag = true) const;

    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const VariableOrderMarkov **imarkov ,
                                                      stat_tool::FrequencyDistribution **hlength ,
                                                      const char *path = NULL) const;
    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const VariableOrderMarkov **markov , int nb_sequence ,
                                                      int length , const char *path = NULL) const;
    stat_tool::DistanceMatrix* divergence_computation(stat_tool::StatError &error , std::ostream &os , int nb_model ,
                                                      const VariableOrderMarkov **markov , int nb_sequence ,
                                                      const MarkovianSequences **seq , const char *path = NULL) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    VariableOrderMarkovData* get_markov_data() const { return markov_data; }
    CategoricalSequenceProcess* get_state_process() const
    { return state_process; }
    int get_nb_output_process() const { return nb_output_process; }
    CategoricalSequenceProcess** get_categorical_process() const
    { return categorical_process; }
    CategoricalSequenceProcess* get_categorical_process(int variable) const
    { return categorical_process[variable]; }
    stat_tool::DiscreteParametricProcess** get_discrete_parametric_process() const
    { return discrete_parametric_process; }
    stat_tool::DiscreteParametricProcess* get_discrete_parametric_process(int variable) const
    { return discrete_parametric_process[variable]; }
    stat_tool::ContinuousParametricProcess** get_continuous_parametric_process() const
    { return continuous_parametric_process; }
    stat_tool::ContinuousParametricProcess* get_continuous_parametric_process(int variable) const
    { return continuous_parametric_process[variable]; }
  };


  VariableOrderMarkov* variable_order_markov_ascii_read(stat_tool::StatError &error ,
                                                        const char *path ,
                                                        int length = DEFAULT_LENGTH);



  class VariableOrderMarkovIterator {  // iterateur chaine de Markov d'ordre variable

  private :

    VariableOrderMarkov *markov;  // pointeur sur un objet VariableOrderMarkov
    int memory;             // memoire

    void copy(const VariableOrderMarkovIterator &it);

  public :

    VariableOrderMarkovIterator(VariableOrderMarkov *imarkov);
    VariableOrderMarkovIterator(const VariableOrderMarkovIterator &iter)
    { copy(iter); }
    ~VariableOrderMarkovIterator();
    VariableOrderMarkovIterator& operator=(const VariableOrderMarkovIterator &iter);

    bool simulation(int **int_seq , int ilength = 1 , bool initialization = false);
    int** simulation(int ilength = 1 , bool initialization = false);

    // acces membres de la classe

    VariableOrderMarkov* get_markov() const { return markov; }
    int get_memory() const { return memory; }
    int get_nb_variable() const { return (markov ? markov->nb_output_process + 1 : 0); }
  };



  class VariableOrderMarkovChainData : public stat_tool::ChainData {  // structure de donnees correspondant
                                                                      // a une chaine de Markov d'ordre variable

  public :

    VariableOrderMarkovChainData(char type , int inb_state , int inb_row , bool init_flag = false)
    :ChainData(type , inb_state , inb_row , init_flag) {}

    void estimation(VariableOrderMarkovChain &markov , bool non_terminal = false ,
                    int estimator = MAXIMUM_LIKELIHOOD ,
                    double laplace_coeff = LAPLACE_COEFF) const;
  };



  class VariableOrderMarkovData : public MarkovianSequences {  // structure de donnees correspondant
                                                               // a une chaine de Markov d'ordre variable

    friend class MarkovianSequences;
    friend class VariableOrderMarkov;
    friend class HiddenVariableOrderMarkov;

    friend std::ostream& operator<<(std::ostream &os , const VariableOrderMarkovData &seq)
    { return seq.ascii_write(os , false); }

  private :

    VariableOrderMarkov *markov;  // pointeur sur un objet VariableOrderMarkov
    VariableOrderMarkovChainData *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences observees
    double restoration_likelihood;  // vraisemblance des sequences restaurees
    double sample_entropy;  // entropie des sequences d'etats
    double *posterior_probability;  // probabilite a posteriori de la sequence d'etats la plus probable
    double *entropy;        // entropie des sequences d'etats
    double *nb_state_sequence;  // nombre de sequences d'etats

    void copy(const VariableOrderMarkovData &seq , bool model_flag = true);
    void observation_frequency_distribution_correction(FrequencyDistribution **corrected_observation ,
                                                       int variable , int start) const;

  public :

    VariableOrderMarkovData();
    VariableOrderMarkovData(const stat_tool::FrequencyDistribution &ihlength , int inb_variable ,
                            int *itype , bool init_flag = false);
    VariableOrderMarkovData(const MarkovianSequences &seq);
    VariableOrderMarkovData(const MarkovianSequences &seq , char transform ,
                            bool initial_run_flag);
    VariableOrderMarkovData(const VariableOrderMarkovData &seq , bool model_flag = true ,
                            char transform = 'c')
    :MarkovianSequences(seq , transform) { copy(seq , model_flag); }
    ~VariableOrderMarkovData();
    VariableOrderMarkovData& operator=(const VariableOrderMarkovData &seq);

    DiscreteDistributionData* extract(stat_tool::StatError &error , int type ,
                                      int variable , int value) const;
    VariableOrderMarkovData* remove_index_parameter(stat_tool::StatError &error) const;
    VariableOrderMarkovData* explicit_index_parameter(stat_tool::StatError &error) const;
    MarkovianSequences* build_auxiliary_variable(stat_tool::StatError &error) const;

    Correlation* state_autocorrelation_computation(stat_tool::StatError &error , int istate ,
                                                   int max_lag = MAX_LAG) const;
    Correlation* output_autocorrelation_computation(stat_tool::StatError &error , int variable ,
                                                    int output , int max_lag = MAX_LAG) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(stat_tool::StatError &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(stat_tool::StatError &error , const char *path , bool exhaustive = false) const;
    bool spreadsheet_write(stat_tool::StatError &error , const char *path) const;
    bool plot_write(stat_tool::StatError &error , const char *prefix , const char *title = NULL) const;
    stat_tool::MultiPlotSet* get_plotable() const;

    void build_transition_count(const VariableOrderMarkovChain &markov ,
                                bool begin = true , bool non_terminal = false);
    void order0_estimation(VariableOrderMarkov &markov) const;

    // acces membres de la classe

    VariableOrderMarkov* get_markov() const { return markov; }
    VariableOrderMarkovChainData* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
    double get_restoration_likelihood() const { return restoration_likelihood; }
    double get_sample_entropy() const { return sample_entropy; }
    double get_posterior_probability(int index) const { return posterior_probability[index]; }
    double get_entropy(int index) const { return entropy[index]; }
    double get_nb_state_sequence(int index) const { return nb_state_sequence[index]; }
  };


};  // namespace sequence_analysis



#endif
