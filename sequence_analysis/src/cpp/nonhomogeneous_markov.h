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
 *       $Id: markov.h 3257 2007-06-06 12:56:12Z dufourko $
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



#ifndef NONHOMOGENEOUS_MARKOV_H
#define NONHOMOGENEOUS_MARKOV_H



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
    friend class Nonhomogeneous_markov;
    friend class Markovian_sequences;
    friend class Nonhomogeneous_markov_data;

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



// class Nonhomogeneous_markov : public STAT_interface , public Chain {
class Nonhomogeneous_markov : public STAT_interface , protected Chain {  // chaine de Markov non-homogene

    friend class Markovian_sequences;
    friend class Nonhomogeneous_markov_data;

    friend Nonhomogeneous_markov* nonhomogeneous_markov_ascii_read(Format_error &error , const char *path ,
                                                                   int length);
    friend std::ostream& operator<<(std::ostream &os , const Nonhomogeneous_markov &markov)
    { return markov.ascii_write(os , markov.markov_data); }

protected :

    Nonhomogeneous_markov_data *markov_data;  // pointeur sur un objet Nonhomogeneous_markov_data
    bool *homogeneity;      // homogeneite des etats
    Function **self_transition;  // fonction d'evolution des probabilites
                                 // de rester dans un etat
    Nonparametric_sequence_process *process;

    void copy(const Nonhomogeneous_markov &markov , bool data_flag = true ,
              bool characteristic_flag = true);
    void remove();

    std::ostream& ascii_write(std::ostream &os , const Nonhomogeneous_markov_data *seq ,
                              bool exhaustive = false , bool file_flag  = false) const;
    std::ostream& spreadsheet_write(std::ostream &os , const Nonhomogeneous_markov_data *seq) const;
    bool plot_write(const char *prefix , const char *title ,
                    const Nonhomogeneous_markov_data *seq) const;

    int nb_parameter_computation() const;

    void transition_update(int state , int index , Chain &index_chain) const;
    void index_state_distribution();
    void state_no_occurrence_probability(int state , double increment = LEAVE_INCREMENT);
    void state_first_occurrence_distribution(int state , int min_nb_value = 1 ,
                                             double cumul_threshold = CUMUL_THRESHOLD);
    void state_nb_pattern_mixture(int state , char pattern);

public :

    Nonhomogeneous_markov();
    Nonhomogeneous_markov(int inb_state , int *ident);
    Nonhomogeneous_markov(const Chain *pchain , const Function **pself_transition , int length);
    Nonhomogeneous_markov(const Nonhomogeneous_markov &markov , bool data_flag = true ,
                          bool characteristic_flag = true)
    :Chain(markov) { copy(markov , data_flag , characteristic_flag); }
    virtual ~Nonhomogeneous_markov();
    Nonhomogeneous_markov& operator=(const Nonhomogeneous_markov &markov);

    Parametric_model* extract(Format_error &error , int type , int state) const;

    std::ostream& line_write(std::ostream &os) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Nonhomogeneous_markov);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void characteristic_computation(int length , bool counting_flag);
    void characteristic_computation(const Nonhomogeneous_markov_data &seq , bool counting_flag ,
                                    bool length_flag = true);

    double likelihood_computation(const Markovian_sequences &seq ,
                                  int index = I_DEFAULT) const;

    Nonhomogeneous_markov_data* simulation(Format_error &error , const Histogram &hlength ,
                                           bool counting_flag = true) const;
    Nonhomogeneous_markov_data* simulation(Format_error &error , int nb_sequence , int length ,
                                           bool counting_flag = true) const;
    Nonhomogeneous_markov_data* simulation(Format_error &error , int nb_sequence ,
                                           const Markovian_sequences &iseq ,
                                           bool counting_flag = true) const;

    // acces membres de la classe

    Nonhomogeneous_markov_data* get_markov_data() const { return markov_data; }
    bool get_homogeneity(int state) const { return homogeneity[state]; }
    Function* get_self_transition(int state) const { return self_transition[state]; }
    Nonparametric_sequence_process* get_process() const { return process; }
};


Nonhomogeneous_markov* nonhomogeneous_markov_ascii_read(Format_error &error , const char *path ,
                                                        int length = DEFAULT_LENGTH);



class Nonhomogeneous_markov_data : public Markovian_sequences {  // structure de donnees correspondant
                                                                 // a une chaine de Markov non-homogene

    friend class Markovian_sequences;
    friend class Nonhomogeneous_markov;

    friend std::ostream& operator<<(std::ostream &os , const Nonhomogeneous_markov_data &seq)
    { return seq.ascii_write(os , false); }

private :

    Nonhomogeneous_markov *markov;  // pointeur sur un objet Nonhomogeneous_markov
    Chain_data *chain_data;  // etats initaux et transitions
    double likelihood;      // vraisemblance des sequences

    void copy(const Nonhomogeneous_markov_data &seq , bool model_flag = true);

public :

    Nonhomogeneous_markov_data();
    Nonhomogeneous_markov_data(const Histogram &ihlength);
    Nonhomogeneous_markov_data(const Markovian_sequences &seq);
    Nonhomogeneous_markov_data(const Nonhomogeneous_markov_data &seq , bool model_flag = true ,
                               char transform = 'c')
    :Markovian_sequences(seq , transform) { copy(seq , model_flag); }
    ~Nonhomogeneous_markov_data();
    Nonhomogeneous_markov_data& operator=(const Nonhomogeneous_markov_data &seq);

    Distribution_data* extract(Format_error &error , int type , int state) const;
    Nonhomogeneous_markov_data* remove_index_parameter(Format_error &error) const;

    std::ostream& ascii_write(std::ostream &os , bool exhaustive = false) const;
    bool ascii_write(Format_error &error , const char *path ,
                     bool exhaustive = false) const;
    bool spreadsheet_write(Format_error &error , const char *path) const;
    bool plot_write(Format_error &error , const char *prefix ,
                    const char *title = 0) const;

/*    RWDECLARE_COLLECTABLE(Nonhomogeneous_markov_data);

    RWspace binaryStoreSize() const;
    void restoreGuts(RWvistream&);
    void restoreGuts(RWFile&);
    void saveGuts(RWvostream&) const;
    void saveGuts(RWFile&) const; */

    void build_transition_count();

    // acces membres de la classe

    Nonhomogeneous_markov* get_markov() const { return markov; }
    Chain_data* get_chain_data() const { return chain_data; }
    double get_likelihood() const { return likelihood; }
};



#endif
