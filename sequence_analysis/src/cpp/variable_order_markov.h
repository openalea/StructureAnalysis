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



#ifndef VARIABLE_ORDER_MARKOV_H
#define VARIABLE_ORDER_MARKOV_H



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


// class Variable_order_markov : public STAT_interface , public Chain {
class Variable_order_markov : public STAT_interface , protected Chain {  // chaine de Markov
                                                                         // d'ordre variable

    friend class Markovian_sequences;
    friend class Variable_order_markov_iterator;
    friend class Variable_order_markov_data;
    friend class Variable_order_chain_data;

    friend Variable_order_markov* variable_order_markov_parsing(Format_error &error ,
                                                                ifstream &in_file ,
                                                                int &line , char type);
    friend Variable_order_markov* variable_order_markov_ascii_read(Format_error &error ,
                                                                   const char *path ,
                                                                   int length);
    friend std::ostream& operator<<(std::ostream &os , const Variable_order_markov &markov)
    { return markov.ascii_write(os); }

protected :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    Variable_order_markov_data *markov_data;  // pointeur sur un objet Variable_order_markov_data
    int *memory_type;       // types des memoires (NON_TERMINAL/TERMINAL/COMPLETION)
    int *order;             // ordres des memoires
    int max_order;          // ordre maximum des memoires
    int **state;            // etats composant les memoires
    int *parent;            // memoires parent
    int **child;            // memoires fils
    int **next;             // memoires suivantes
    int *nb_memory;         // nombre de memoires precedentes
    int **previous;         // memoires precedentes
    int nb_output_process;  // nombre de processus d'observation
    Nonparametric_sequence_process **nonparametric_process;  // processus d'observation non-parametriques
    Parametric_process **parametric_process;  // processus d'observation parametriques

    Variable_order_markov(const Variable_order_markov *pmarkov , int inb_output_process ,
                          Nonparametric_process **nonparametric_observation ,
                          Parametric_process **parametric_observation , int length);

    void memory_tree_completion(const Variable_order_markov &markov);
    void copy(const Variable_order_markov &markov , bool data_flag = true);
    void remove();

    std::ostream& ascii_memory_tree_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& ascii_transition_tree_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , const Variable_order_markov_data *seq ,
                              bool exhaustive = false , bool file_flag = false ,
                              bool hidden = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Variable_order_markov_data *seq ,
                                    bool hidden = false) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Variable_order_markov_data *seq) const;

    void find_parent_memory(int index);
    void build_memory_transition();
    void build_previous_memory();
    bool check_free_suffix() const;
    bool** logic_transition_computation() const;
    void component_computation();

    void build_non_terminal();

    void threshold_application(double min_probability);

    void max_order_computation();
    int nb_parameter_computation(double min_probability = 0.) const;
    int nb_transient_parameter_computation(double min_probability = 0.) const;
    double penalty_computation(bool hidden , double min_probability = 0.) const;

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

    Correlation* state_autocorrelation_computation(Format_error &error ,
                                                   int istate , int max_lag ,
                                                   const Variable_order_markov_data *seq) const;
    Correlation* output_autocorrelation_computation(Format_error &error , int variable ,
                                                    int output , int max_lag ,
                                                    const Variable_order_markov_data *seq) const;

    double likelihood_computation(const Variable_order_chain_data &chain_data) const;

    double likelihood_correction(const Variable_order_markov_data &seq) const;

    std::ostream& transition_count_ascii_write(std::ostream &os , bool begin) const;

public :

    Variable_order_markov();
    Variable_order_markov(char itype , int inb_state , int inb_row);
    Variable_order_markov(char itype , int inb_state , int inb_row ,
                          int imax_order);
    Variable_order_markov(char itype , int inb_state , int iorder , bool init_flag ,
                          int inb_output_process = 0 , int nb_value = 0);
    Variable_order_markov(const Variable_order_markov &markov ,
                          int inb_output_process , int nb_value);
/*    Variable_order_markov(const Variable_order_markov &markov ,
                          int inb_output_process , int *nb_value); */
    Variable_order_markov(const Variable_order_markov *pmarkov ,
                          const Nonparametric_process *pobservation , int length);
    Variable_order_markov(const Variable_order_markov &markov , bool data_flag = true)
    :Chain(markov) { copy(markov , data_flag); }
    virtual ~Variable_order_markov();
    Variable_order_markov& operator=(const Variable_order_markov &markov);

    Parametric_model* extract(Format_error &error , int type ,
                              int variable , int value) const;
    Variable_order_markov_data* extract_data(Format_error &error) const;

    Variable_order_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    // fonctions pour la compatibilite avec la classe STAT_interface

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const Variable_order_markov_data &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    Correlation* state_autocorrelation_computation(Format_error &error , int istate ,
                                                   int max_lag = MAX_LAG) const;
    Correlation* output_autocorrelation_computation(Format_error &error , int variable ,
                                                    int output , int max_lag = MAX_LAG) const;

    double likelihood_computation(const Markovian_sequences &seq , int index) const;
    double likelihood_computation(const Variable_order_markov_data &seq) const;

    Variable_order_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                           bool counting_flag = true , bool divergence_flag = false) const;
    Variable_order_markov_data* simulation(Format_error &error , int nb_sequence ,
                                           int length , bool counting_flag = true) const;
    Variable_order_markov_data* simulation(Format_error &error , int nb_sequence ,
                                           const Markovian_sequences &iseq ,
                                           bool counting_flag = true) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Variable_order_markov **imarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Variable_order_markov **markov , int nb_sequence ,
                                            int length , const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Variable_order_markov **markov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    Variable_order_markov_data* get_markov_data() const { return markov_data; }
    int get_memory_type(int memory) const { return memory_type[memory]; }
    int get_order(int memory) const { return order[memory]; }
    int get_max_order() const { return max_order; }
    int get_state(int memory , int lag) const { return state[memory][lag]; }
    int get_parent(int memory) const { return parent[memory]; }
    int get_child(int memory , int istate) const { return child[memory][istate]; }
    int get_next(int memory , int istate) const { return next[memory][istate]; }
    int get_nb_memory(int memory) const { return nb_memory[memory]; }
    int get_previous(int memory , int istate) const { return previous[memory][istate]; }
    int get_nb_output_process() const { return nb_output_process; }
    Nonparametric_sequence_process* get_nonparametric_process(int variable) const
    { return nonparametric_process[variable]; }
    Parametric_process** get_parametric_process() const { return parametric_process; }
    Parametric_process* get_parametric_process(int variable)
    const { return parametric_process[variable]; }
};


Variable_order_markov* variable_order_markov_parsing(Format_error &error ,
                                                     ifstream &in_file ,
                                                     int &line , char type);
Variable_order_markov* variable_order_markov_ascii_read(Format_error &error ,
                                                        const char *path ,
                                                        int length = DEFAULT_LENGTH);



class Variable_order_markov_iterator {  // iterateur chaine de Markov d'ordre variable

private :

    Variable_order_markov *markov;  // pointeur sur un objet Variable_order_markov
    int memory;             // memoire

    void copy(const Variable_order_markov_iterator &it);

public :

    Variable_order_markov_iterator(Variable_order_markov *imarkov);
    Variable_order_markov_iterator(const Variable_order_markov_iterator &it)
    { copy(it); }
    ~Variable_order_markov_iterator();
    Variable_order_markov_iterator& operator=(const Variable_order_markov_iterator &it);

    bool simulation(int **int_seq , int ilength = 1 , bool initialization = false);
    int** simulation(int ilength = 1 , bool initialization = false);

    // acces membres de la classe

    Variable_order_markov* get_markov() const { return markov; }
    int get_memory() const { return memory; }
    int get_nb_variable() const { return (markov ? markov->nb_output_process + 1 : 0); }
};



class Variable_order_chain_data : public Chain_data {  // structure de donnees correspondant
                                                       // a une chaine de Markov d'ordre variable

public :

    Variable_order_chain_data(char type , int inb_state , int inb_row , bool init_flag = false)
    :Chain_data(type , inb_state , inb_row , init_flag) {}

    void estimation(Variable_order_markov &markov , bool non_terminal = false ,
                    int estimator = MAXIMUM_LIKELIHOOD ,
                    double laplace_coeff = LAPLACE_COEFF) const;
};



class Variable_order_markov_data : public Markovian_sequences {  // structure de donnees correspondant
                                                                 // a une chaine de Markov d'ordre variable

    friend class Markovian_sequences;
    friend class Variable_order_markov;
    friend class Hidden_variable_order_markov;

    friend std::ostream& operator<<(std::ostream &os , const Variable_order_markov_data &seq)
    { return seq.ascii_write(os , false); }

private :

    Variable_order_markov *markov;  // pointeur sur un objet Variable_order_markov
    Variable_order_chain_data *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences
    double hidden_likelihood;  // vraisemblance de toutes les sequences possibles
    double *posterior_probability;  // probabilite a posteriori de la sequence d'etats la plus probable

    void copy(const Variable_order_markov_data &seq , bool model_flag = true);
    void observation_histogram_correction(Histogram **corrected_observation ,
                                          int variable , int start) const;

public :

    Variable_order_markov_data();
    Variable_order_markov_data(const Histogram &ihlength , int inb_variable , bool init_flag = false);
    Variable_order_markov_data(const Markovian_sequences &seq);
    Variable_order_markov_data(const Markovian_sequences &seq , char transform ,
                               bool initial_run_flag);
    Variable_order_markov_data(const Variable_order_markov_data &seq , bool model_flag = true ,
                               char transform = 'c')
    :Markovian_sequences(seq , transform) { copy(seq , model_flag); }
    ~Variable_order_markov_data();
    Variable_order_markov_data& operator=(const Variable_order_markov_data &seq);

    Distribution_data* extract(Format_error &error , int type ,
                               int variable , int value) const;
    Variable_order_markov_data* remove_index_parameter(Format_error &error) const;

    Correlation* state_autocorrelation_computation(Format_error &error , int istate ,
                                                   int max_lag = MAX_LAG) const;
    Correlation* output_autocorrelation_computation(Format_error &error , int variable ,
                                                    int output , int max_lag = MAX_LAG) const;

    std::ostream& ascii_data_write(std::ostream &os , char format = 'c' ,
                                   bool exhaustive = false) const;
    bool ascii_data_write(Format_error &error , const char *path ,
                          char format = 'c' , bool exhaustive = false) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

    void build_transition_count(const Variable_order_markov &markov ,
                                bool begin = true , bool non_terminal = false);
    void order0_estimation(Variable_order_markov &markov) const;

    // acces membres de la classe

    Variable_order_markov* get_markov() const { return markov; }
    Variable_order_chain_data* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
    double get_hidden_likelihood() const { return hidden_likelihood; }
    double get_posterior_probability(int index) const { return posterior_probability[index]; }
};



#endif
