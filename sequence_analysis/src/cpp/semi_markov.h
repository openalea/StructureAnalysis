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


// class Semi_markov : public STAT_interface , public Chain {
class Semi_markov : public STAT_interface , protected Chain {  // semi-chaine de Markov

    friend class Markovian_sequences;
    friend class Semi_markov_iterator;
    friend class Semi_markov_data;

    friend Semi_markov* semi_markov_ascii_read(Format_error &error , const char *path ,
                                               int length, bool counting_flag ,
                                               double cumul_threshold );
    friend std::ostream& operator<<(std::ostream &os , const Semi_markov &smarkov)
    { return smarkov.ascii_write(os , smarkov.semi_markov_data); }

protected :

    int nb_iterator;        // nombre d'iterateurs pointant sur l'objet
    Semi_markov_data *semi_markov_data;  // pointeur sur un objet Semi_markov_data
    int *state_subtype;     //  MARKOVIAN/SEMI_MARKOVIAN
    Forward **forward;      // lois de l'intervalle de temps residuel
    int nb_output_process;  // nombre de processus d'observation
    Nonparametric_sequence_process **nonparametric_process;  // processus d'observation non-parametriques
    Parametric_process **parametric_process;  // processus d'observation parametriques

    Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                int inb_output_process , Nonparametric_process **pobservation ,
                int length , bool counting_flag);
    Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                int inb_output_process , Nonparametric_process **nonparametric_observation ,
                Parametric_process **parametric_observation , int length , bool counting_flag);

    void copy(const Semi_markov &smarkov , bool data_flag = true ,
              int param = I_DEFAULT);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Semi_markov_data *seq ,
                              bool exhaustive = false , bool file_flag  = false ,
                              bool hidden = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Semi_markov_data *seq ,
                                    bool hidden = false) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Semi_markov_data *seq) const;

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

    Semi_markov();
    Semi_markov(char itype , int inb_state , int inb_output_process , int *nb_value);
    Semi_markov(const Chain *pchain , const Nonparametric_sequence_process *poccupancy ,
                const Nonparametric_process *pobservation , int length ,
                bool counting_flag);
    Semi_markov(const Semi_markov &smarkov , bool data_flag = true ,
                int param = I_DEFAULT)
    :Chain(smarkov) { copy(smarkov , data_flag , param); }
    void conditional_delete();
    virtual ~Semi_markov();
    Semi_markov& operator=(const Semi_markov &smarkov);

    Parametric_model* extract(Format_error &error , int type ,
                              int variable , int value) const;
    Parametric_model* extract(Format_error &error , int state ,
                              int histogram_type = FINAL_RUN) const;
    Semi_markov_data* extract_data(Format_error &error) const;

    Semi_markov* thresholding(double min_probability = MIN_PROBABILITY) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Semi_markov);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void characteristic_computation(int length , bool counting_flag , int variable = I_DEFAULT);
    void characteristic_computation(const Semi_markov_data &seq , bool counting_flag ,
                                    int variable = I_DEFAULT , bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq , int index) const;
    double likelihood_computation(const Semi_markov_data &seq) const;

    Semi_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                 bool counting_flag = true , bool divergence_flag = false) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence , int length ,
                                 bool counting_flag = true) const;
    Semi_markov_data* simulation(Format_error &error , int nb_sequence ,
                                 const Markovian_sequences &iseq , bool counting_flag = true) const;

    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **ismarkov , Histogram **hlength ,
                                            const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **smarkov , int nb_sequence ,
                                            int length , const char *path = 0) const;
    Distance_matrix* divergence_computation(Format_error &error , std::ostream &os , int nb_model ,
                                            const Semi_markov **smarkov , int nb_sequence ,
                                            const Markovian_sequences **seq , const char *path = 0) const;

    // acces membres de la classe

    int get_nb_iterator() const { return nb_iterator; }
    Semi_markov_data* get_semi_markov_data() const { return semi_markov_data; }
    int get_state_subtype(int state) const { return state_subtype[state]; }
    Forward** get_forward() const { return forward; }
    Forward* get_forward(int state) const { return forward[state]; }
    int get_nb_output_process() const { return nb_output_process; }
    Nonparametric_sequence_process* get_nonparametric_process(int variable)
    const { return nonparametric_process[variable]; }
    Parametric_process** get_parametric_process() const { return parametric_process; }
    Parametric_process* get_parametric_process(int variable)
    const { return parametric_process[variable]; }
};


Semi_markov* semi_markov_ascii_read(Format_error &error , const char *path ,
                                    int length = DEFAULT_LENGTH , bool counting_flag = true ,
                                    double cumul_threshold = OCCUPANCY_THRESHOLD);



class Semi_markov_iterator {  // iterateur semi-chaine de Markov

private :

    Semi_markov *semi_markov;  // pointeur sur un objet Semi_markov
    int state;              // etat
    int occupancy;          // temps d'occupation de l'etat
    int counter;            // compteur

    void copy(const Semi_markov_iterator &it);

public :

    Semi_markov_iterator(Semi_markov *ismarkov);
    Semi_markov_iterator(const Semi_markov_iterator &it)
    { copy(it); }
    ~Semi_markov_iterator();
    Semi_markov_iterator& operator=(const Semi_markov_iterator &it);

    void simulation(int **seq , int length = 1 , bool initialization = false);
    int** simulation(int length = 1 , bool initialization = false);

    // acces membres de la classe

    Semi_markov* get_semi_markov() const { return semi_markov; }
    int get_state() const { return state; }
    int get_occupancy() const { return occupancy; }
    int get_counter() const { return counter; }
    int get_nb_variable() const { return (semi_markov ? semi_markov->nb_output_process + 1 : 0); }
};



class Semi_markov_data : public Markovian_sequences {  // structure de donnees correspondant
                                                       // a une semi-chaine de Markov

    friend class Markovian_sequences;
    friend class Semi_markov;
    friend class Hidden_semi_markov;

    friend std::ostream& operator<<(std::ostream &os , const Semi_markov_data &seq)
    { return seq.ascii_write(os , false); }

private :

    Semi_markov *semi_markov;  // pointeur sur un objet Semi_markov
    Chain_data *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences
    double hidden_likelihood;  // vraisemblance de toutes les sequences possibles
    double *posterior_probability;  // probabilite a posteriori de la sequence d'etats la plus probable

    void copy(const Semi_markov_data &seq , bool model_flag = true);

public :

    Semi_markov_data();
    Semi_markov_data(const Histogram &ihlength , int inb_variable , bool init_flag = false);
    Semi_markov_data(const Markovian_sequences &seq);
    Semi_markov_data(const Markovian_sequences &seq , char transform , bool initial_run_flag);
    Semi_markov_data(const Semi_markov_data &seq , bool model_flag = true , char transform = 'c')
    :Markovian_sequences(seq , transform) { copy(seq , model_flag); }
    ~Semi_markov_data();
    Semi_markov_data& operator=(const Semi_markov_data &seq);

    Distribution_data* extract(Format_error &error , int type ,
                               int variable , int value) const;
    Semi_markov_data* remove_index_parameter(Format_error &error) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Semi_markov_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void build_transition_count(const Semi_markov *smarkov = 0);

    // acces membres de la classe

    Semi_markov* get_semi_markov() const { return semi_markov; }
    Chain_data* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
    double get_hidden_likelihood() const { return hidden_likelihood; }
    double get_posterior_probability(int index) const { return posterior_probability[index]; }
};



#endif
