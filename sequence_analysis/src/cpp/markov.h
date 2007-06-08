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



#ifndef MARKOV_H
#define MARKOV_H



/****************************************************************
 *
 *  Constantes :
 */


const double START_RATIO = 0.03;       // proportion de l'echantillon pour l'initialisation
                                       // des parametres (debut)
const double END_RATIO = 0.1;          // proportion de l'echantillon pour l'initialisation
                                       // des parametres (fin)
const int REGRESSION_NB_ELEMENT = 100;  // taille minimum de l'echantillon pour
                                        // la regression non-lineare
const double GRADIENT_DESCENT_COEFF = 1.;  // coefficient algorithme de gradient
const double RESIDUAL_SQUARE_SUM_DIFF = 1.e-6;  // seuil pour stopper les iterations
                                                // de l'algorithme de gradient
const int REGRESSION_NB_ITER = 1000;   // nombre d'iterations pour la regression non-lineaire



/****************************************************************
 *
 *  Definition des classes :
 */


class Function : public Regression_kernel {  // fonction d'evolution des probabilites
                                             // de rester dans un etat

    friend class Self_transition;
    friend class Markov;
    friend class Markov_iterator;
    friend class Markovian_sequences;
    friend class Markov_data;

    friend Function* function_parsing(Format_error &error , std::ifstream &in_file , int &line ,
                                      int length , double min, double max);

private :

    double *residual;       // residus
    int *frequency;         // effectifs correspondant a chaque index

    void copy(const Function&);
    void remove();

    std::ostream& ascii_print(std::ostream &os , bool exhaustive , bool file_flag ,
                              const Curves *curves = 0) const;
    std::ostream& spreadsheet_print(std::ostream &os , const Curves *curves = 0) const;
    bool plot_print(const char *path , double residual_standard_deviation = D_DEFAULT) const;

    double regression_square_sum_computation(double self_transition_mean) const;
    void residual_computation(const Self_transition &self_transition);
    double residual_mean_computation() const;
    double residual_variance_computation(double residual_mean) const;
    double residual_square_sum_computation() const;

public :

    Function();
    Function(int iident , int length , double *iparameter);
    Function(int iident , int length);
    Function(const Function &function);
    ~Function();
    Function& operator=(const Function &function);

/*    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    // acces membres de la classe

    double get_residual(int index) const { return residual[index]; }
    int get_frequency(int index) const { return frequency[index]; }
};


Function* function_parsing(Format_error &error , std::ifstream &in_file , int &line ,
                           int length , double min = 0. , double max = 1.);



// class Markov : public STAT_interface , public Chain {
class Markov : public STAT_interface , protected Chain {  // chaine de Markov

    friend class Markovian_sequences;
    friend class Markov_iterator;
    friend class Markov_data;

    friend Markov* markov_ascii_read(Format_error &error , const char *path ,
                                     int length);
    friend std::ostream& operator<<(std::ostream &os , const Markov &markov)
    { return markov.ascii_write(os , markov.markov_data); }

protected :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    Markov_data *markov_data;  // pointeur sur un objet Markov_data
    int order;              // ordre
    int *self_row;          // indices des lignes de la matrice des probabilites de transition
                            // correspondant aux probabilites de rester dans un etat
    bool homogeneous;       // homogene/non-homogene
    bool *homogeneity;      // homogeneite des etats
    Function **self_transition;  // fonction d'evolution des probabilites
                                 // de rester dans un etat
    int nb_output_process;  // nombre de processus d'observation
    Nonparametric_sequence_process **process;  // processus d'observation non-parametriques

    Markov(const Chain *pchain , int iorder , int inb_output_process ,
           Nonparametric_process **pobservation , int length);
    Markov(const Markov &markov , int state);

    void copy(const Markov &markov , bool data_flag = true ,
              bool characteristic_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Markov_data *seq ,
                              bool exhaustive = false , bool file_flag  = false ,
                              bool hidden = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Markov_data *seq ,
                                    bool hidden = false) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Markov_data *seq) const;

    void self_row_computation();
    void component_computation();

    int nb_parameter_computation(double min_probability = 0.) const;

    void index_state_distribution();
    double* memory_computation() const;
    void state_no_occurrence_probability(int state , double increment = LEAVE_INCREMENT);
    void state_first_occurrence_distribution(int state , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void state_leave_probability(const double *memory , int state ,
                                 double increment = LEAVE_INCREMENT);
    void state_recurrence_time_distribution(const double *memory , int state ,
                                            int min_nb_value = 1 ,
                                            double cumul_threshold = CUMUL_THRESHOLD);
    void state_sojourn_time_distribution(const double *memory , int state ,
                                         int min_nb_value = 1 ,
                                         double cumul_threshold = CUMUL_THRESHOLD);
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

    double likelihood_correction(const Markov_data &seq) const;

    // chaine de Markov non-homogene

    void transition_update(int state , int index , Chain &index_chain) const;
    void non_homogeneous_index_state_distribution();
    void non_homogeneous_state_no_occurrence_probability(int state ,
                                                         double increment = LEAVE_INCREMENT);
    void non_homogeneous_state_first_occurrence_distribution(int state , int min_nb_value = 1 ,
                                                             double cumul_threshold = CUMUL_THRESHOLD);
    void non_homogeneous_state_nb_pattern_mixture(int state , char pattern);

public :

    Markov();
    Markov(int inb_state , int *ident);
    Markov(int inb_state , int iorder , int inb_output_process , int *nb_value);
    Markov(const Chain *pchain , int iorder , const Function **pself_transition ,
           const Nonparametric_process *pobservation , int length);
    Markov(int iorder , const Markov &markov);
    Markov(const Markov &markov , bool data_flag = true ,
           bool characteristic_flag = true)
    :Chain(markov) { copy(markov , data_flag , characteristic_flag); }
    virtual ~Markov();
    Markov& operator=(const Markov &markov);

    Parametric_model* extract(Format_error &error , int type ,
                              int variable , int value) const;
    Markov_data* extract_data(Format_error &error) const;

    Markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    Markov* MTD_model_building(Format_error &error , std::ostream &os , int iorder ,
                               double *weight , int length = DEFAULT_LENGTH) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Markov);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const Markov_data &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq , int index) const;
    double likelihood_computation(const Markov_data &seq) const;

    Markov_data* simulation(Format_error &error , const Histogram &hlength ,
                            bool counting_flag = true , bool divergence_flag = false) const;
    Markov_data* simulation(Format_error &error , int nb_sequence , int length ,
                            bool counting_flag = true) const;
    Markov_data* simulation(Format_error &error , int nb_sequence ,
                            const Markovian_sequences &iseq , bool counting_flag = true) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Markov **imarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Markov **markov , int nb_sequence , int length ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Markov **markov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;

    // chaine de Markov non-homogene

    void non_homogeneous_characteristic_computation(int length , bool counting_flag);
    void non_homogeneous_characteristic_computation(const Markov_data &seq , bool counting_flag ,
                                                    bool length_flag = true);

    double non_homogeneous_likelihood_computation(const Markovian_sequences &seq ,
                                                  int index = I_DEFAULT) const;

    Markov_data* non_homogeneous_simulation(Format_error &error , const Histogram &hlength ,
                                            bool counting_flag = true) const;
    Markov_data* non_homogeneous_simulation(Format_error &error , int nb_sequence , int length ,
                                            bool counting_flag = true) const;
    Markov_data* non_homogeneous_simulation(Format_error &error , int nb_sequence ,
                                            const Markovian_sequences &iseq ,
                                            bool counting_flag = true) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    Markov_data* get_markov_data() const { return markov_data; }
    int get_order() const { return order; }
    int get_self_row(int state) const { return self_row[state]; }
    bool get_homogeneous() const { return homogeneous; }
    bool get_homogeneity(int state) const { return homogeneity[state]; }
    Function* get_self_transition(int state) const { return self_transition[state]; }
    int get_nb_output_process() const { return nb_output_process; }
    Nonparametric_sequence_process* get_process(int variable) const
    { return process[variable]; }
};


Markov* markov_ascii_read(Format_error &error , const char *path ,
                          int length = DEFAULT_LENGTH);



class Markov_iterator {  // iterateur chaine de Markov

private :

    Markov *markov;         // pointeur sur un objet Markov
    Chain *chain;           // chaine de Markov courante
    int *state;             // etats precedents

    void copy(const Markov_iterator &it);

public :

    Markov_iterator(Markov *imarkov);
    Markov_iterator(const Markov_iterator &it)
    { copy(it); }
    ~Markov_iterator();
    Markov_iterator& operator=(const Markov_iterator &it);

    void simulation(int **seq , int length = 1 , bool initialization = false);
    int** simulation(int length = 1 , bool initialization = false);

    // acces membres de la classe

    Markov* get_markov() const { return markov; }
    Chain* get_chain() const { return chain; }
    int get_state(int index) const { return state[index]; }
    int get_nb_variable() const { return (markov ? markov->nb_output_process + 1 : 0); }
};



class Markov_data : public Markovian_sequences {  // structure de donnees correspondant
                                                  // a une chaine de Markov

    friend class Markovian_sequences;
    friend class Markov;
    friend class Hidden_markov;

    friend std::ostream& operator<<(std::ostream &os , const Markov_data &seq)
    { return seq.ascii_write(os , false); }

private :

    Markov *markov;         // pointeur sur un objet Markov
    Chain_data *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences
    double hidden_likelihood;  // vraisemblance de toutes les sequences possibles

    void copy(const Markov_data &seq , bool model_flag = true);

public :

    Markov_data();
    Markov_data(int inb_variable , const Histogram &ihlength , bool init_flag);
    Markov_data(int inb_sequence , int *ilength , int ***isequence ,
                const Markov &imarkov);
    Markov_data(const Markovian_sequences &seq);
    Markov_data(const Markovian_sequences &seq , int variable);
    Markov_data(const Markov_data &seq , bool model_flag = true)
    :Markovian_sequences(seq) { copy(seq , model_flag); }
    ~Markov_data();
    Markov_data& operator=(const Markov_data &seq);

    Distribution_data* extract(Format_error &error , int type ,
                               int variable , int value) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Markov_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void build_transition_count(int order , bool begin = true);

    // acces membres de la classe

    Markov* get_markov() const { return markov; }
    Chain_data* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
    double get_hidden_likelihood() const { return hidden_likelihood; }
};



#endif
