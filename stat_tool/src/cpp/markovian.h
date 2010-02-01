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



#ifndef MARKOVIAN_H
#define MARKOVIAN_H



#include "chain_reestimation.h"



/****************************************************************
 *
 *  Constantes :
 */


const int NB_STATE = 100;              // nombre maximum d'etats d'une chaine de Markov
const int ORDER = 8;                   // ordre maximum d'une chaine de Markov
const double MIN_PROBABILITY = 1.e-5;  // probabilite minimum
const double THRESHOLDING_FACTOR = 0.8;  // facteur pour le seuillage des probabilites
const int NB_PARAMETER = 100000;       // nombre maximum de parametres d'une chaine de Markov
const int NB_OUTPUT_PROCESS = 10;      // nombre maximum de processus d'observation
const int NB_OUTPUT = 12;              // nombre maximum d'observations par etat (cas non-parametrique)
const double OBSERVATION_THRESHOLD = 0.999;  // seuil sur la fonction de repartition
                                             // pour borner une loi d'observation parametrique

const double ACCESSIBILITY_THRESHOLD = 1.e-6;  // seuil pour stopper l'algorithme
                                               // de calcul de l'accessibilite des etats
const int ACCESSIBILITY_LENGTH = 100;  // longueur maximum de sequence pour l'algorithme
                                       // de calcul de l'accessibilite des etats

const double NOISE_PROBABILITY = 0.05;  // perturbation des probabilites d'observation

const int MIN_NB_ELEMENT = 10;         // taille minimum de l'echantillon construit par arrondi
const int OBSERVATION_COEFF = 10;      // coefficient arrondi estimateur pour les lois
                                       // d'observation parametriques

// const double SELF_TRANSITION = 0.9;    probabilite de rester dans un etat initiale

enum {
  FORWARD ,
  FORWARD_BACKWARD ,
  VITERBI ,
//  VITERBI_FORWARD_BACKWARD ,
  GENERALIZED_VITERBI ,
  FORWARD_BACKWARD_SAMPLING ,
  GIBBS_SAMPLING ,
  FORWARD_DYNAMIC_PROGRAMMING
};



/****************************************************************
 *
 *  Definition des classes :
 */


class Chain_data;

class Chain {           // chaine de Markov

//    friend class Chain_data;

    friend Chain* chain_parsing(Format_error &error , ifstream &in_file ,
                                int &line , char type);
    friend std::ostream& operator<<(std::ostream &os , const Chain &chain)
    { return chain.ascii_print(os , true); }

// protected :
public :

    char type;              // 'o' : ordinaire, 'e' : en equilibre
    int nb_state;           // nombre d'etats
    int nb_row;             // nombre de lignes de la matrice
                            // des probabilites de transition
    bool **accessibility;   // matrice d'accessibilite des etats
    int nb_component;       // nombre de classes
    int *component_nb_state;  // nombre d'etats par classe
    int **component;        // classes
    char *state_type;       // types des etats ('r' : recurrent,
                            // 't' : transitoire, 'a' : absorbant)
    double *initial;        // probabilites initiales
    double *cumul_initial;  // fonction de repartition correspondant
                            // au probabilites initiales
    double **transition;    // matrice des probabilites de transition
    double **cumul_transition;  // fonctions de repartition correspondant aux lignes
                                // de la matrice des probabilites de transition

    void parameter_copy(const Chain&);
    void copy(const Chain&);
    void remove();

    std::ostream& ascii_print(std::ostream &os , bool file_flag = false) const;
    std::ostream& spreadsheet_print(std::ostream &os) const;

    void create_cumul();
    void cumul_computation();
    void remove_cumul();
    void log_computation();

    bool** logic_transition_computation() const;
    bool connex_component_research(Format_error &error , bool **ilogic_transition = NULL) const;
    void graph_accessibility_computation(bool **ilogic_transition);
    void probability_accessibility_computation();
    void component_computation(bool **ilogic_transition = NULL);

    void thresholding(double min_probability);

    int nb_parameter_computation(double min_probability = 0.) const;
    double chi2_value_computation(const Chain_data &chain_data) const;

// public :

    Chain(char itype = 'o' , int inb_state = 0 , bool init_flag = true);
    Chain(char itype , int inb_state , int inb_row , bool init_flag);
    Chain(const Chain &chain) { copy(chain); }
    ~Chain();
    Chain& operator=(const Chain &chain);

    void init(bool left_right , double self_transition);

    double likelihood_computation(const Chain_data &chain_data , bool initial_flag = true) const;
    void chi2_fit(const Chain_data &chain_data , Test &test) const;

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_nb_state() const { return nb_state; }
    int get_nb_row() const { return nb_row; }
    bool get_accessibility(int state1 , int state2) const
    { return accessibility[state1][state2]; }
    int get_nb_component() const { return nb_component; }
    int get_component_nb_state(int icomponent) const
    { return component_nb_state[icomponent]; }
    int get_component(int icomponent , int index) const
    { return component[icomponent][index]; }
    char get_state_type(int state) const { return state_type[state]; }
    double get_initial(int state) const { return initial[state]; }
    double get_transition(int memory , int state) const
    { return transition[memory][state]; } */
};


Chain* chain_parsing(Format_error &error , ifstream &in_file , int &line , char type);



class Chain_data : public Chain_reestimation<int> {  // structure de donnees correspondant a
                                                     // une chaine de Markov

/*    friend class Chain;
    friend class Markovian_sequences;
    friend class Variable_order_markov;
    friend class Variable_order_markov_data;
    friend class Semi_markov_data;
    friend class Non_homogeneous_markov;
    friend class Non_homogeneous_markov_data; */

// protected :
public :

    int nb_parameter_computation() const;

// public :

    Chain_data(char itype , int inb_state , int inb_row , bool init_flag = false)
    :Chain_reestimation<int>(itype , inb_state , inb_row , init_flag) {}
    Chain_data(const Chain_data &chain_data);

    void estimation(Chain &chain) const;

    // acces membres de la classe

/*    char get_type() const { return type; }
    int get_nb_state() const { return nb_state; }
    int get_nb_row() const { return nb_row; }
    int get_initial(int state) const { return initial[state]; }
    int get_transition(int memory , int state) const
    { return transition[memory][state]; } */
};



class Nonparametric_process {  // processus d'observation non-parametrique

    friend class Mv_Mixture;

    friend Nonparametric_process* observation_parsing(Format_error &error , ifstream &in_file ,
                                                      int &line , int nb_state , bool hidden);

    friend Nonparametric_process** observation_parsing(Format_error &error , ifstream &in_file ,
                                                       int &line , int nb_state , int &nb_output_process);
    friend Nonparametric_process** old_observation_parsing(Format_error &error , ifstream &in_file ,
                                                           int &line , int nb_state , int &nb_output_process);

protected :

    int nb_state;           // nombre d'etats
    int nb_value;           // nombre de valeurs
    Distribution **observation;  // lois d'observation

    void copy(const Nonparametric_process &process);
    void add_state(const Nonparametric_process &process , int state);
    void remove();

    bool test_hidden() const;
    void thresholding(double min_probability);

    void state_permutation(int *perm) const; // permutation des etats

public :

    Nonparametric_process(int inb_state = 0 , int inb_value = 0 , int observation_flag = false);
    Nonparametric_process(int inb_state , int inb_value , double **observation_probability);
    Nonparametric_process(int inb_state , Distribution **pobservation);
    Nonparametric_process(const Nonparametric_process &process , char manip = 'c' , int state = I_DEFAULT);
    ~Nonparametric_process();
    Nonparametric_process& operator=(const Nonparametric_process &process);

    int nb_parameter_computation(double min_probability) const;
    void init();

    // acces membres de la classe

    int get_nb_state() const { return nb_state; }
    int get_nb_value() const { return nb_value; }
    Distribution* get_observation(int state) const
    { return observation[state]; }

    std::ostream& ascii_print(std::ostream &os , Histogram **empirical_observation ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os ,
                                    Histogram **empirical_observation = NULL) const;
    bool plot_print(const char *prefix , const char *title ,
                    int process , Histogram **empirical_observation = NULL) const;
};


Nonparametric_process* observation_parsing(Format_error &error , ifstream &in_file ,
                                           int &line , int nb_state , bool hidden);

Nonparametric_process** observation_parsing(Format_error &error , ifstream &in_file ,
                                            int &line , int nb_state , int &nb_output_process);
Nonparametric_process** old_observation_parsing(Format_error &error , ifstream &in_file ,
                                                int &line , int nb_state , int &nb_output_process);



class Parametric_process {  // processus d'observation parametrique

/*    friend class Markovian_sequences;
    friend class Semi_markov;
    friend class Semi_markov_data;
    friend class Semi_markov_iterator;
    friend class Hidden_semi_markov;
    friend class Variable_order_markov;
    friend class Variable_order_markov_data;
    friend class Variable_order_markov_iterator;
    friend class Hidden_variable_order_markov; */

    friend Parametric_process* observation_parsing(Format_error &error , ifstream &in_file ,
                                                   int &line , int nb_state ,
                                                   double cumul_threshold);

// protected :
public :

    int nb_state;           // nombre d'etats
    int nb_value;           // nombre de valeurs
    Parametric **observation;  // lois d'observation

    void copy(const Parametric_process &process);

    void state_permutation(int *perm) const; // permutation des etats

    void add_state(const Parametric_process &process , int state);
    void remove();

    void nb_value_computation();

// public :

    Parametric_process(int inb_state = 0 , int inb_value = 0);
    Parametric_process(int inb_state , Parametric **pobservation);
    Parametric_process(const Parametric_process &process , char manip = 'c' , int state = I_DEFAULT);
    ~Parametric_process();
    Parametric_process& operator=(const Parametric_process &process);

    std::ostream& ascii_print(std::ostream &os , Histogram **empirical_observation ,
                              bool exhaustive , bool file_flag) const;
    std::ostream& spreadsheet_print(std::ostream &os , Histogram **empirical_observation = NULL) const;
    bool plot_print(const char *prefix , const char *title , int process ,
                    Histogram **empirical_observation = NULL) const;
    void plotable_write(MultiPlotSet &plot , int &index , int process ,
                        Histogram **empirical_observation = NULL) const;

    int nb_parameter_computation() const;
    void init();

    // acces membres de la classe

/*    int get_nb_state() const { return nb_state; }
    int get_nb_value() const { return nb_value; }
    Parametric* get_observation(int state) const
    { return observation[state]; } */
};


Parametric_process* observation_parsing(Format_error &error , ifstream &in_file ,
                                        int &line , int nb_state ,
                                        double cumul_threshold = OBSERVATION_THRESHOLD);



# endif
